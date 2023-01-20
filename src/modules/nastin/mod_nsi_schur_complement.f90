!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_schur_complement.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   This routine solves iterations of the Schur complement
!> @details This routine solves iterations of the Schur complement
!>          solver Richardson - Orthomin(1)
!-----------------------------------------------------------------------

module mod_nsi_schur_complement

  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_nastin
  use mod_postpr
  use def_kermod,              only : kfl_duatss
  use mod_maths,               only : maths_equalize_arrays
  use mod_memory,              only : memory_alloca_min
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_solver,              only : solver_parallel_SpMV
  use mod_local_basis,         only : local_basis_global_to_local
  use mod_local_basis,         only : local_basis_local_to_global
  use mod_nsi_schur_operations
  use mod_couplings,           only : couplings_check_dirichlet 
  use mod_messages,            only : livinf
  use mod_communications,      only : PAR_BARRIER
  use mod_communications,      only : PAR_SUM
  use mod_arrays,              only : arrays_number
  use mod_run_config,          only : run_config
  implicit none

  private
  
  integer(ip), private                   :: idofn,ipoin,nu,np,nzP
  integer(ip), private                   :: ipass
  real(rp),    private                   :: numer,denom
  real(rp),    private                   :: ak=1.0_rp,raux2,rauxc,rauxm
  real(rp),                      pointer :: xx(:)
  real(rp),                      pointer :: zz(:)
  real(rp),                      pointer :: vv(:)
  real(rp),                      pointer :: rr(:)
  real(rp),                      pointer :: bp(:)
  real(rp),                      pointer :: bu(:)
  real(rp),                      pointer :: uu(:)
  real(rp),                      pointer :: pp(:)
  real(rp),                      pointer :: Q(:)
  real(rp),                      pointer :: zzold(:)
  real(rp),                      pointer :: Tim(:)     ! Nodal projection of dt/rho
  real(rp),                      pointer :: Tau(:)     ! Nodal projection of tau
#ifdef outmateocoe
  integer(ip)        :: ipass_aux
#endif
  real(rp),parameter                     :: gama = 1.0_rp  ! leave this always =1.0 (putting 0.0 we can recover some non incremental versions)

  public :: nsi_schur_complement_initialization
  public :: nsi_schur_complement_solution
  public :: nsi_schur_complement_memory
  public :: nsi_schur_complement_matrices
  public :: nsi_momentum_continuity_residuals
  public :: xx
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Initialization
  !> @details Initialization of Schur complement method
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_schur_complement_initialization()

    ipass = 0
    nullify(xx)
    nullify(zz)
    nullify(vv)
    nullify(rr)
    nullify(bp)
    nullify(bu)
    nullify(uu)
    nullify(pp)
    nullify(Q)
    nullify(zzold)
    nullify(Tim)   
    nullify(Tau)   

  end subroutine nsi_schur_complement_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Assembly and solution
  !> @details Assembly of the NS equations and solution using the
  !>          different solution method for the pressure Schur complement
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_schur_complement_solution()

    if ( kfl_expco_nsi /= 1 .or. itinn(modul) == 1 ) then  ! for the case kfl_expco_nsi==1 do this part only the first iteration, else do it always

       !-------------------------------------------------------------------
       !                 +-       -+ +- -+   +-  -+
       !                 | Auu Aup | | u |   | bu |
       ! Assemble system |         | |   | = |    |
       !                 | Apu App | | p |   | bp |
       !                 +-       -+ +- -+   +-  -+
       !-------------------------------------------------------------------

       call nsi_solini(3_ip)                                        ! Initialize solver
       call local_basis_global_to_local(kfl_fixrs_nsi,uu)           ! Global to local
       call nsi_matrix()                                            ! Assemble equation

       !-------------------------------------------------------------------
       !                                            ----------------
       ! Schur complement preconditioning Q = App - Apu Auu^{-1} Aup
       ! if it was not done previously
       !
       !-------------------------------------------------------------------

       call nsi_solini(2_ip)
       call nsi_schpre()

    end if

    !-------------------------------------------------------------------
    !
    ! Iterative Schur complement solver
    !
    !-------------------------------------------------------------------

    call livinf(56_ip,' ',modul)    
    call livinf(160_ip,' ',1_ip)

    if( NSI_SCHUR_COMPLEMENT ) then

       if( kfl_sosch_nsi == 2 ) then                 ! Orthomin(1): momentum preserving
          call nsi_solmom()
          call nsi_solort()
       else if( kfl_sosch_nsi == 3 ) then            ! Orthomin(1): continuity preserving
          call nsi_solmom()
          call nsi_solort()
          call nsi_solric()
       else if( kfl_sosch_nsi == 4 ) then            ! Richardson
          call nsi_solmom()
          call nsi_solric()
       else if( kfl_sosch_nsi == 5 ) then            ! Richardson: momentum preserving
          call nsi_solmom()
          call nsi_solric()
          call nsi_solmom()
       else if( kfl_sosch_nsi == 6 ) then            ! Richardson: continuity preserving
          call nsi_solmom()
          call nsi_solric()
       end if

       call local_basis_local_to_global(kfl_fixrs_nsi,uu)          ! Local to global
       !
       ! Velocity correction
       !
       if( kfl_sosch_nsi == 3 .or. kfl_sosch_nsi == 6 ) then
          if (kfl_corre_nsi == 3) then
             call nsi_solcor_consistent_mass()
          else
             call nsi_solcor()
          end if
       end if

    end if

    call livinf(164_ip,' ',1_ip)

  end subroutine nsi_schur_complement_solution

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Solve momentum eqns
  !> @details Momentum equation solution for Richardson and Orthomin(1)
  !>          methods
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_solmom()

    !-------------------------------------------------------------------
    !
    ! Momentum equation
    !
    ! INPUT  p^k
    ! OUTPUT u^k, r^k
    !
    !-------------------------------------------------------------------
    implicit none

    call times(1) % ini() 
    !
    ! Solve the momentum equation: Auu u^{k+1} = bu - Aup p^k
    !
    nmome_nsi = nmome_nsi + 1
    call nsi_solini(1_ip)                                                ! Initialize solver

    if( INOTEMPTY ) then
       if (nint(gama) == 1) then
          call nsi_aupvec(1_ip,amatr(poaup_nsi:),unkno(ndbgs_nsi+1:),vv)  ! vv = Aup p^k
       else
          vv = 0.0_rp
       end if
       do idofn = 1,nu                                                   ! vv = bu - Aup p^k
          vv(idofn) = bu(idofn) - vv(idofn)
       end do

    end if

    call times(1) % add() 

    if( kfl_duatss == 1 ) then
       call nsi_sol_dts_prec(vv,uu,amatr(poauu_nsi),pmatr)
    else
       !
       ! For the moment I will only leave the second velocity solved with AGMG
       ! For the second velocity the initial guess es 0 and therefore ratio and toler
       ! are actually the same thing. We put 2 different values but actually the one
       ! that counts is the smaller (typically ratio). AGMG only has stopping by toler.
       ! Therefore Alya's solver and AGMG ae only comparable if we do not use adaptive
       ! (which would significantly alter our typical startegy) or if the initial guess is 0.
       !
!#ifdef solve_w_agmg
!       call nsi_agmgsol(vv,uu,amatr(poauu_nsi),r_sol,c_sol,ndime,npoin)
!#else
       call solver(vv,uu,amatr(poauu_nsi),pmatr)                          ! Solve system
!#endif
  
    end if

    call livinf(165_ip,'M',0_ip)

   
    call times(1) % ini() 
    if( INOTEMPTY ) then

       call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,rr)                   ! rr =  Apu u^{k+1}

       do ipoin = 1,np                                                 ! rr <= bp - rr
          rr(ipoin) = rhsid(ndbgs_nsi+ipoin) - rr(ipoin)
       end do
       if (nint(gama) == 1) call nsi_appvec(2_ip,amatr(poapp_nsi:),unkno(ndbgs_nsi+1:),rr)  ! rr <= rr - App p^k

    end if
 
    call times(1) % add() 

  end subroutine nsi_solmom

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Schur complement system with a preconditioned Richardson
  !>          iteration
  !> @details Solve and actualize pressure
  !>          INPUT  r^k
  !>          OUTPUT p^{k+1}
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_solric()

    call times(1) % ini() 
    
    nschu_nsi = nschu_nsi + 1
    call nsi_solini(2_ip)                                     ! Initialize solver
    !
    ! Change residual norm
    !
    if( solve_sol(1) % kfl_normalization == 1 ) then
       call nsi_appvec(0_ip,Q,pp,vv)
       call rhsmod(1_ip,vv)
       call norm2x(1_ip,vv,solve_sol(1) % normalization)
    end if
    !
    ! Solve: Q z = r^k   with   Q = App + P
    !
    call nsi_inivec(np,zz)                                    ! z := Dp = 0
    ! For gama==0 I could use  a better initialization
    !
    ! Solve system
    !
    call times(1) % add() 
    call solver(rr,zz,Q,pmatr)                                ! Solve system

    call livinf(165_ip,'S',0_ip)
    call times(1) % ini() 
    !
    ! Actualize p^{k+1}
    !
    ! Q Dp' = rc
    !
    !
    ! p^{k+1} = Dp'+p^k
    !
    do ipoin = 1,np
       unkno(ndbgs_nsi+ipoin) = zz(ipoin) + unkno(ndbgs_nsi+ipoin)
    end do

    if (nint(gama) == 1) then  ! Normal Alya behaviour
       do ipoin = 1,np
          unkno(ndbgs_nsi+ipoin) = zz(ipoin) + unkno(ndbgs_nsi+ipoin)
       end do
    else   ! notice that in nsi_elmcor no modification is necesary because in deltp_nsi =>zz it will hav P not deltaP
       do ipoin = 1,np
          unkno(ndbgs_nsi+ipoin) = zz(ipoin)
       end do
    end if

    call times(1) % add() 
    
  end subroutine nsi_solric

  subroutine nsi_solcor()

    !-------------------------------------------------------------------
    !
    ! Correct velocity
    !
    ! INPUT  u^{k}, p^{k+1}, p^k
    ! OUTPUT u^{k+1}
    !
    !-------------------------------------------------------------------

    !                                  ----------
    ! Correction M u^{k+1} = M u^{k} - Auu^-1 Aup [p^{k+1})-p^k]
    !
    call nsi_elmcor()
    if(run_config%timing) call PAR_BARRIER()

    call livinf(165_ip,'C',0_ip)

  end subroutine nsi_solcor

  !-----------------------------------------------------------------------
  !> 
  !> @author  Herbert Owen
  !> @date    2018-05-10
  !> @brief   Correct velocity using consistent mass. Moreover the correction uses Aup directly
  !>          not an approximation Auu^-1 Aup by the Galerkin part of Aup.
  !> @details See Algebraic pressure segregation methods for the incompressible Navier-Stokes equations S. Badia · R. Codina
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_solcor_consistent_mass()
    !
    ! Correction M ( u^{k+1} - u^{k} )  = - Aup [p^{k+1})-p^k]
    !
    !
    ! From nsi_solric we have delta p in zz
    ! Now do xx = - Aup zz - Since Aup is in local coordinates xx is also in local coordinates
    !
    call nsi_aupvec(1_ip,amatr(poaup_nsi:),zz,xx)       ! xx = Aup zz
    !
    ! Beware since it is defined as xx(:)  I can not do xx = -xx -
    ! the same applies latter for uu = uu + vv
    ! I believe this is NOT NICE
    !
    xx(1:ndime*npoin) = -xx(1:ndime*npoin)

    call nsi_inisol(5_ip)
    call solver(xx,vv,cmama_nsi,pmatr)                 ! Solve system  I put delta_u in vv
    ! Here I leave pmatr as is done in the rest of the subroutine but actually it is not used so you need not care.
    ! It would be used if you put PRECONDITIONER: MATRIX  (see reasol and all_matrix). We typically use DIAGO or LINEL. This I believe is not nice.
    ! What would happen if somebody put PRECONDITIONER: MATRIX in the nsi.dat file??

    call local_basis_local_to_global(kfl_fixrs_nsi,vv)           ! Local to global for delta_u

    uu(1:ndime*npoin) = uu(1:ndime*npoin) + vv(1:ndime*npoin)

    if(run_config%timing) call PAR_BARRIER()    ! par_barrier

    call livinf(165_ip,'C',0_ip)

  end subroutine nsi_solcor_consistent_mass

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Schur complement system with a preconditioned Orthomin(1) iteration
  !> @details Schur complement system with a preconditioned Orthomin(1) iteration
  !>          INPUT  u^{k+1}
  !>          OUTPUT u^{k+2}, p^{k+1}, r^{k+1}
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_solort()
    !
    ! Solve: Q z = r^k   with   Q = App + P
    !
    nmome_nsi = nmome_nsi + 1
    nschu_nsi = nschu_nsi + 1
    call nsi_solini(2_ip)                                     ! Initialize solver
    call times(1) % ini() 
    !
    ! Change residual norm
    !
    if( solve_sol(1) % kfl_normalization == 1 ) then
       call nsi_appvec(0_ip,Q,pp,vv)
       call rhsmod(1_ip,vv)
       call norm2x(1_ip,vv,solve_sol(1) % normalization)
    end if
 
    if( kfl_incre_nsi == 0 ) then
       call solver_parallel_SpMV(solve(2),lapla_nsi,unkno(ndbgs_nsi+1:),zz,MPI=.false.,OPENMP=.true.)
       do ipoin = 1,npoin
          vv(ipoin) = rr(ipoin)
          rr(ipoin) = rr(ipoin) + zz(ipoin)
          zz(ipoin) = unkno(ndbgs_nsi+ipoin)
       end do
    else
       call nsi_inivec(npoin,zz)                              ! Dp=0
    end if
    !
    ! Solve system
    !
#ifdef outmateocoe
    call nsi_matndof(rr,zz,lapla_nsi,r_dom,c_dom,1_ip,npoin,ipass_aux)
#endif

    call times(1) % add() 
    if( mitri_nsi == 1 ) then
#ifdef solve_w_agmg
       call nsi_agmgsol(rr,zz,lapla_nsi,r_dom,c_dom,1_ip,npoin)
#else
       call solver(rr,zz,lapla_nsi,pmatr)                        ! Solve system
#endif
    else
       call nsi_richar()
    end if
      
#ifdef outmateocoe
    call nsi_out_res_eocoe(zz,npoin,ipass_aux)
#endif

    call times(1) % ini() 

    call livinf(165_ip,'S',0_ip)

    if( kfl_incre_nsi == 0 ) then
       do ipoin = 1,npoin
          zz(ipoin) = zz(ipoin) - unkno(ndbgs_nsi+ipoin)
          rr(ipoin) = vv(ipoin)
       end do
    end if
    !
    ! Solve: Auu v = Aup z
    !
    call nsi_solini(4_ip)
    call nsi_aupvec(1_ip,amatr(poaup_nsi:),zz,xx)             ! xx = Aup z
    call nsi_inivec(nu,vv)
    
    call times(1) % add() 

#ifdef outmateocoe
    !
    ! transform matrix to agmg format and output
    !
    call nsi_matndof(xx,vv,amatr(poauu_nsi),r_sol,c_sol,ndime,npoin,ipass_aux)
#endif

    if( kfl_duatss == 1 ) then
       call nsi_sol_dts_prec(xx,vv,amatr(poauu_nsi),pmatr)
    else
#ifdef solve_w_agmg
       call nsi_agmgsol(xx,vv,amatr(poauu_nsi),r_sol,c_sol,ndime,npoin)
#else
       call solver(xx,vv,amatr(poauu_nsi),pmatr)               ! Solve system
#endif
    end if
        

#ifdef outmateocoe
    call nsi_out_res_eocoe(vv,ndime*npoin,ipass_aux)
#endif

    call times(1) % ini() 
    call livinf(165_ip,'M',0_ip)
    !
    ! Compute x = App z - Apu v
    !
    call nsi_appvec(1_ip,amatr(poapp_nsi:),zz,xx)              ! xx = App z
    call nsi_apuvec(2_ip,amatr(poapu_nsi:),vv,xx)              ! rr = Apu v
    call rhsmod(1_ip,xx)
    !
    ! Compute ak = <rr,xx>/<xx,xx>
    !
    !call solver_parallel_scalar_product(solve(2),xx,xx,rr,denom,numer)
    call prodts(1_ip,np,xx,rr,denom,numer) 
    if( denom == 0.0_rp ) then
       ak = 1.0_rp
    else
       ak = numer / denom
    end if
    !
    ! Actualize p^{k+1} = p^k     + ak*z
    !           u^{k+2} = u^{k+1} - ak*v
    !
    do ipoin = 1,np
       unkno(ndbgs_nsi+ipoin) = unkno(ndbgs_nsi+ipoin) + ak * zz(ipoin)
    end do
    
    do idofn = 1,nu
       uu(idofn) = uu(idofn) - ak * vv(idofn)
    end do
    !
    ! Actualize r^{k+1} = r^k - ak*xx
    !
    call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,rr)                  ! rr = Apu u^{k+1}
    do ipoin = 1,np                                                ! rr = bp - rr
       rr(ipoin) = rhsid(ndbgs_nsi+ipoin) - rr(ipoin)
    end do

    call nsi_appvec(2_ip,amatr(poapp_nsi:),unkno(ndbgs_nsi+1:),rr) ! rr = rr - App p^k

    call times(1) % add() 

  end subroutine nsi_solort
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Compute or point to matrices 
  !> @details Compute or point to matrices and vectors
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_schur_complement_matrices()

    use mod_memory
    
   if( IMASTER ) then                                         ! Pointers
       nu  =  0
       np  =  0
       bu  => nul1r
       bp  => nul1r
       uu  => nul1r
       pp  => nul1r
       Q   => nul1r
       nzp =  0
    else
       nu  =  ndime*npoin
       np  =  npoin
       bu  => rhsid
       bp  => rhsid(ndbgs_nsi+1:)
       uu  => unkno
       pp  => unkno(ndbgs_nsi+1:)
       Q   => lapla_nsi
       nzp =  size(Q,KIND=ip)
    end if

    Tim       => dt_rho_nsi
    Tau       => tau_nsi
    deltp_nsi => zz

    !print*,'a=',kfl_paral,memory_size(deltp_nsi),memory_size(Q),memory_size(bu)

  end subroutine nsi_schur_complement_matrices
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Initialize pointers 
  !> @details Initialize pointers 
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_schur_complement_memory()
     
    nu  = ndime*npoin
    np  = npoin
    
    call memory_alloca(mem_modul(1:2,modul),'XX','mod_nsi_solsch',xx,max(1_ip,nu))
    call memory_alloca(mem_modul(1:2,modul),'ZZ','mod_nsi_solsch',zz,max(1_ip,nu))
    call memory_alloca(mem_modul(1:2,modul),'VV','mod_nsi_solsch',vv,max(1_ip,nu))
    call memory_alloca(mem_modul(1:2,modul),'RR','mod_nsi_solsch',rr,max(1_ip,np))    
 
  end subroutine nsi_schur_complement_memory
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   Schur complement and momentum residuals
  !> @details Compute momentum and continuity residuals for
  !>          convergence and postprocess purposes:
  !>          Continuity: rp = bp - Apu.u - App.p
  !>          Momentum:   ru = bu - Auu.u - Aup.p
  !>
  !>          OUTPUT resin_nsi, resou_nsi or resss_nsi
  !>
  !----------------------------------------------------------------------

  subroutine nsi_momentum_continuity_residuals(itask)
    integer(ip), intent(in) :: itask
    integer(ip)             :: ii
    real(rp)                :: resim,resis

    if( itask == 1 ) then                                               ! Inner residual
       ii = 1
    else if( itask == 2 ) then                                          ! Outer residual
       ii = 2
       if( momod(modul) % miinn == 1 .and. micou(iblok) <= 1 ) then
          resou_nsi(1) = resin_nsi(1)
          resou_nsi(2) = resin_nsi(2)
          return
       end if
    else                                                                ! Steady state residual
       ii = nprev_nsi
       if( momod(modul) % miinn == 1 .and. micou(iblok) <= 1 ) then
          resss_nsi(1) = resin_nsi(1)
          resss_nsi(2) = resin_nsi(2)
          return
       end if
    end if
    !
    ! Gobal to local
    !
    call local_basis_global_to_local(kfl_fixrs_nsi,veloc,LAST_COMPONENT=ii)                         ! Global to local
    !
    ! Momentum residual: || [bu - Aup p^k] - Auu u^k || // || bu - Aup p^k ||
    !
    call nsi_aupvec(1_ip,amatr(poaup_nsi:),press(:,ii),vv)          ! vv = Aup p^k
    do ipoin = 1,nu                                                 ! vv = bu - vv
       xx(ipoin) = bu(ipoin) - vv(ipoin)
       vv(ipoin) = xx(ipoin)
    end do
    call nsi_auuvec(2_ip,amatr(poauu_nsi:),veloc(:,:,ii),xx)       ! vv = vv - Auu u^k

    call nsi_solini(1_ip)
    call rhsmod(ndime,vv)
    call rhsmod(ndime,xx)
    call norm2x(ndime,vv,rauxm)
    call norm2x(ndime,xx,raux2)

    if( rauxm <= zeror ) rauxm = 1.0_rp
    resim = raux2 / (rauxm+zeror)                                   ! Residual norm for u
#ifdef DETAILS_ORTHOMIN
    if(kfl_paral<=1) write(777,'(a,3(e9.2))')'VEL:resim,raux2,rauxm',resim,raux2,rauxm
#endif
    if( postp(1) % npp_stepi(arrays_number('MOMEN'),1) > 0 .and. INOTMASTER ) then
       call maths_equalize_arrays(xx,remom_nsi)
    end if
    !
    ! Continuity residual: || [bp - Apu u^k] - App p^k || // || bp - Apu u^k ||
    !
    call nsi_apuvec(1_ip,amatr(poapu_nsi:),veloc(:,:,ii),rr)        ! rr = Apu u^k
    do ipoin = 1,np                                                 ! rr = bp - rr
       zz(ipoin) = bp(ipoin) - rr(ipoin)
       rr(ipoin) = zz(ipoin)
    end do
    call nsi_appvec(2_ip,amatr(poapp_nsi:),press(:,ii),zz)          ! rr = rr - App p^k

    call nsi_solini(2_ip)
    call rhsmod(1_ip,rr)
    call rhsmod(1_ip,zz)
    call norm2x(1_ip,rr,rauxc)
    call norm2x(1_ip,zz,raux2)

    if( rauxc <= zeror ) rauxc = 1.0_rp
    resis = raux2 / (rauxc+zeror)                                   ! Residual norm for p
#ifdef DETAILS_ORTHOMIN
    if(kfl_paral<=1) write(777,'(a,3(e9.2))'),'PRES:resis,raux2,rauxc',resis,raux2,rauxc
#endif
    if( postp(1) % npp_stepi(arrays_number('SCHUR'),1) > 0  .and. INOTMASTER ) then
       call maths_equalize_arrays(zz,resch_nsi)
    end if
    !
    ! Local to global
    !
    call local_basis_local_to_global(kfl_fixrs_nsi,veloc,LAST_COMPONENT=ii)                         ! Local to global

    if(itask==1) then
       resin_nsi(1)=resim
       resin_nsi(2)=resis
    else if(itask==2) then
       resou_nsi(1)=resim
       resou_nsi(2)=resis
    else
       resss_nsi(1)=resim
       resss_nsi(2)=resis
    end if

  end subroutine nsi_momentum_continuity_residuals

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Richardson method
  !> @details Solve Q z = r using a Richardson method, where Q=D is the discrete Laplacian
  !>
  !>          z^{i+1} = z^i + C^-1 ( r - D z^i )
  !>
  !>          1. Compute: r2 = r - Q z^i
  !>          2. Solve:   C y = r2 => y = D ^-1 r2
  !>          3. Update:  z^{i+1} = z^i + y
  !>
  !>          Typically, C is the continuous Laplacian (symmetric) and D the discrete one
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_richar()

    use mod_matrix
    implicit none
    integer(ip)          :: iiter
    real(rp),    pointer :: yy(:)
    real(rp),    pointer :: zz_new(:)
    real(rp),    pointer :: r2(:)
    real(rp)             :: resid(2),alpha
    real(rp),    pointer :: invAuu(:)

    if( INOTMASTER ) then
       !
       ! Allocate memory
       !
       allocate( yy(npoin) )
       allocate( zz_new(npoin) )
       allocate( r2(npoin) )
       allocate( invAuu(ndime*npoin) )
       !
       ! Diagonal of Auu
       !
       call matrix_invdia(&
            npoin,ndime,solve(2) % kfl_symme,r_dom,c_dom,&
            amatr(poauu_nsi:),invAuu,memit)

       do ipoin = 1,npoin
          zz_new(ipoin) = zz(ipoin)
          yy(ipoin)     = 0.0_rp
          r2(ipoin)     = 0.0_rp
       end do

    else

       allocate( yy(1) )
       allocate( r2(1) )

    end if

    do iiter = 1,mitri_nsi
       !
       ! [ App - Apu diag(Auu)^-1 Aup ] z^i
       !
       call nsi_shuvec(&
            0_ip,amatr(poauu_nsi),amatr(poapp_nsi),amatr(poapu_nsi),invAuu,amatr(poaup_nsi),&
            zz_new,r2)
       if( INOTMASTER ) then
          do ipoin = 1,npoin
             r2(ipoin) = rr(ipoin) - r2(ipoin)
             yy(ipoin) = 0.0_rp
          end do
       end if
       !
       ! Solve C y = r2
       !
       call solver(r2,yy,lapla_nsi,pmatr)

       alpha = 1.0_rp
       !call nsi_shuvec(&
       !     0_ip,amatr(poapp_nsi),amatr(poapu_nsi),invAuu,amatr(poaup_nsi),&
       !     yy,zz)
       !call prodxy(1_ip,npoin,zz,r2,alpha)
       !call norm2x(1_ip,zz,denom)
       !alpha = alpha / (denom*denom)
       if( INOTMASTER ) then
          resid(1) = 0.0_rp
          resid(2) = 0.0_rp
          do ipoin = 1,npoin
             resid(1)      = resid(1)      + yy(ipoin) * yy(ipoin)
             resid(2)      = resid(2)      + zz_new(ipoin) * zz_new(ipoin)
             zz_new(ipoin) = zz_new(ipoin) + alpha * yy(ipoin)
          end do
       end if

       call PAR_SUM(2_ip,resid)
       resid(1) = sqrt(resid(1)/(resid(2)+zeror))

    end do
    !
    ! Continuity residual MUST be global
    !
    call pararr('SLX',NPOIN_TYPE,npoin,rr)
    !
    ! Update solution and deallocate memory
    !
    if( INOTMASTER ) then
       do ipoin = 1,npoin
          zz(ipoin) = zz_new(ipoin)
       end do
       deallocate( invAuu )
       deallocate( zz_new )
    end if
    deallocate( yy )
    deallocate( r2 )

  end subroutine nsi_richar

end module mod_nsi_schur_complement
!> @}
