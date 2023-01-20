!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_sol_dts_prec.f90
!> @date    17/04/2015
!> @author  Herbert Owen
!> @brief   Replaces the call to solver for the momentum equation by a precondtioned Richardson iteration 
!> @details The preconditioner used is the same matrix but with the transient term increased fact_duatss times.
!!          For the increse a Lumped mass matrix is used.
!!          The Richardson iteration does a maximum of 2*fact_duatss iterartions.  
!!          A^{*}  delta_u^{m} = -Au^{m}+b
!!          with:  A^{*} = ( (factduatss-1) * (1/delta_t) M_L ) + (1/delta_t) (M) + K
!!           ( (factduatss-1) * (1/delta_t) M_L )  has been built by nsi_elmmat  and assembled into lumma
!!          It stops when |-Ax^{m}+b| reaches the tolerance imposed to the solver or it has been reduced 
!!          by a ratio if ADAPTIVE is used. That is identical to what is done for a typical solver.
!!          Perhaps the difference (to the safe side) is that here it is evaluated on teh non preconditioned residuals.
!!          For the internal solvers the same tol and ratio are used.
!> @} 
!------------------------------------------------------------------------
subroutine nsi_sol_dts_prec(rhsid,unkno,amatr,pmatr)
  use def_kintyp,         only :  ip,rp
  use def_master,         only :  kfl_paral, INOTMASTER,modul,IPARALL,nparr,parre
  use def_master,         only :  mem_modul,solve,lumma
  use def_domain,         only :  npoin,nzdom,ndime,r_dom,c_dom
  use def_kermod,         only :  fact_duatss

  use mod_memory,         only :  memory_alloca,memory_alloca_min
  use mod_memory,         only :  memory_deallo
  use mod_iofile,         only :  iofile
  use mod_messages, only : messages_live
  use mod_nsi_schur_operations
  implicit none
  real(rp),    intent(inout)   :: unkno(npoin*ndime)
  real(rp),    intent(inout)   :: amatr(ndime,ndime,nzdom)
  real(rp),    intent(in)      :: pmatr(*)          ! It is not used
  real(rp),    intent(inout)   :: rhsid(npoin*ndime)

  real(rp),pointer             :: amatr_dts(:,:,:)
  real(rp),pointer             :: delta_dts(:)  ! delta_x^{m}
  real(rp),pointer             :: resi(:),resi_auxi(:)      ! residual
  !  real(rp),pointer            :: rhs_dts(:)   ! no need ew can just use resi

  real(rp)                     :: raux0,raux,max_delta
  character(200)               :: auxch

  real(rp),    target          :: rmaxi(1)
  integer(ip)                  :: ipoin,iiter,jpoin,izdom

  nullify(amatr_dts)
  nullify(delta_dts)
  nullify(resi)
  nullify(resi_auxi)
  !     nullify(rhs_dts)

  if( INOTMASTER ) then
     call memory_alloca(mem_modul(1:2,modul),'AMATR_DTS','nsi_sol_dts_prec',amatr_dts,ndime,ndime,nzdom)
     call memory_alloca(mem_modul(1:2,modul),'DELTA_DTS','nsi_sol_dts_prec',delta_dts,npoin*ndime)
     call memory_alloca(mem_modul(1:2,modul),'RESI','nsi_sol_dts_prec',resi,npoin*ndime)
     call memory_alloca(mem_modul(1:2,modul),'RESI_AUXI','nsi_sol_dts_prec',resi_auxi,npoin*ndime)
     !        call memory_alloca(mem_modul(1:2,modul),'RHS_DTS','nsi_sol_dts_prec',rhs_dts,npoin*ndime)
  else
     call memory_alloca_min(mem_modul(1:2,modul),'AMATR_DTS','nsi_sol_dts_prec',amatr_dts)
     call memory_alloca_min(mem_modul(1:2,modul),'DELTA_DTS','nsi_sol_dts_prec',delta_dts)
     call memory_alloca_min(mem_modul(1:2,modul),'RESI'     ,'nsi_sol_dts_prec',resi)
     call memory_alloca_min(mem_modul(1:2,modul),'RESI_AUXI','nsi_sol_dts_prec',resi_auxi)
     !        call memory_alloca_min(rhs_dts)
  end if
  if( INOTMASTER ) then
     !
     ! Obtain initial residual
     !
     call nsi_auuvec(1_ip,amatr,unkno,resi)            ! resi0 = A u^{m=0}
     do ipoin = 1,npoin*ndime                                                 
        resi(ipoin) = rhsid(ipoin) - resi(ipoin)       ! resi0 = rhsid - A u^{m=0}
     end do

     resi_auxi = resi
     call rhsmod(ndime,resi_auxi)
     call norm2x(ndime,resi_auxi,raux0)     

     if(raux0 < solve(1) % solco) then                ! tolerance for momentum
        if (kfl_paral==2) write(*,*) 'initial residual less than tolerance,raux,tol',raux0,solve(1) % solco
     end if
     !
     ! Obtain the modified matrix - A^{*} = ( (factduatss-1) * (1/delta_t) M_L ) + (1/delta_t) (M) + K
     !
     amatr_dts = amatr
     do ipoin = 1,npoin
        izdom = r_dom(ipoin) - 1
        jpoin = 0
        do while( jpoin /= ipoin )
           izdom = izdom + 1
           jpoin = c_dom(izdom)
        end do
        amatr_dts(1,1,izdom) = amatr_dts(1,1,izdom) + lumma(ipoin)
        amatr_dts(2,2,izdom) = amatr_dts(2,2,izdom) + lumma(ipoin)
        if (ndime == 3) amatr_dts(3,3,izdom) = amatr_dts(3,3,izdom) + lumma(ipoin)
     end do
     
  else
     resi_auxi = 0.0_rp
     call norm2x(ndime,resi_auxi,raux0)
  end if

  main_do: do iiter=1,100*fact_duatss   ! ojo este valor de 100 es algo con lo que estoy jugando supongo habría que entrar
     delta_dts = 0.0_rp          ! initialization = 0 because we are solving for the increment
     call solver(resi,delta_dts,amatr_dts,pmatr)                          ! Solve system

!     print*,'kfl_paral,ipoin,lninv_loc(ipoin),delta_dts'

     if( INOTMASTER ) then
        unkno = unkno + delta_dts
        !
        ! Obtain residual 
        !        
        call nsi_auuvec(1_ip,amatr,unkno,resi)            ! resi = A u^m
        do ipoin = 1,npoin*ndime                                                 
           resi(ipoin) = rhsid(ipoin) - resi(ipoin)       ! resi = rhsid - A u^m
        end do

        
        resi_auxi = resi
        call rhsmod(ndime,resi_auxi)
        call norm2x(ndime,resi_auxi,raux)
        max_delta = maxval(abs(delta_dts))
     else  ! to have something for the master - not sure if it is needed
        max_delta = 0.0_rp
        resi_auxi = 0.0_rp
        call norm2x(ndime,resi_auxi,raux)
     end if  ! NOTMASTER

     
  
     if( IPARALL ) then
        nparr     =  1
        rmaxi(1)  =  max_delta
        parre     => rmaxi
        call par_operat(2_ip)
        max_delta =  rmaxi(1)
     end if

     write(auxch,'((a),i5,3(1x,e12.5))') 'iiter,raux,raux/raux0,max_delta',iiter,raux,raux/raux0,max_delta
     call messages_live(auxch) 

     if(raux < solve(1) % solco) exit main_do
     if(raux/raux0 < solve(1) % adres) exit main_do
  end do main_do

  call memory_deallo(mem_modul(1:2,modul),'AMATR_DTS','nsi_sol_dts_prec',amatr_dts)
  call memory_deallo(mem_modul(1:2,modul),'DELTA_DTS','nsi_sol_dts_prec',delta_dts)
  call memory_deallo(mem_modul(1:2,modul),'RESI','nsi_sol_dts_prec',resi)

end subroutine nsi_sol_dts_prec
