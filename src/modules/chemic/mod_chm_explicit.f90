!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_explicit

!  use def_parame
  use def_kintyp,     only : ip, rp
  use def_master,     only : INOTMASTER, ITASK_BEGINN, ITASK_ENDINN, mem_modul, modul, dtinv, rhsid, unkno, itinn
  use def_domain,     only : npoin, vmass
  use def_chemic,     only : dt_rho_chm, comax_chm, comin_chm, dt_chm, iclaf_chm, iclai_chm, ittot_chm, kfl_goit2_chm,&
                             kfl_spray_chm, kfl_tiacc_chm, nclas_chm, rtpts_chm, iclas_chm, kfl_fixno_chm, bvess_chm
  use mod_solver,     only : solver_lumped_mass_system
  use mod_memory,     only : memory_alloca

  implicit none

  real(rp),   pointer, save :: bt(:)    ! Energy RHS
  real(rp),   pointer, save :: tt(:)    ! Temperature
  real(rp),   pointer, save :: rt(:,:,:)  ! Energy residual
  real(rp),   pointer, save :: Mass(:)  ! lumped mass
  real(rp),   pointer, save :: Tim(:,:,:) ! projection dt / (rho*cp)
  real(rp),            save :: dt(3)

  private

  public :: chm_explicit_solution

contains

   subroutine chm_explicit_allocate()

      integer(ip), save :: ipass = 0

      if( ipass == 0 ) then

         ipass = 1
         nullify(rt)
         nullify(Tim)

         bt       => rhsid
         tt       => unkno
         Mass     => vmass

         if( INOTMASTER ) then
               call memory_alloca(mem_modul(1:2,modul),'RT' ,'chm_explicit_allocate',rt, nclas_chm,npoin,4_ip)
               call memory_alloca(mem_modul(1:2,modul),'TIM','chm_explicit_allocate',Tim,nclas_chm,npoin,4_ip)
         else
            call memory_alloca(mem_modul(1:2,modul),'RT' ,'chm_explicit_allocate',rt, 1_ip,1_ip,4_ip)
            call memory_alloca(mem_modul(1:2,modul),'TIM','chm_explicit_allocate',Tim,1_ip,1_ip,4_ip)
         end if

         rt  = 0.0_rp
         Tim = 0.0_rp

         dt = 1.0_rp

      end if

   end subroutine chm_explicit_allocate


   subroutine chm_explicit_solution()
     integer(ip) :: ipoin,time_order,kpoin
     real(rp)    :: reslo
     real(rp)    :: dt01,dt012,dt12,alpha1, alpha2, alpha3
     real(rp)    :: rholoc
     integer(ip), save :: iter = 0

     external    :: chm_updunk
     external    :: chm_matrix
     external    :: chm_endite
     external    :: rhsmod

     iter = iter + 1

     !
     ! Update inner iteration counter
     !
     itinn(modul) = itinn(modul) + 1_ip
     ittot_chm    = ittot_chm + 1_ip
     rtpts_chm     =  0.0_rp
     comin_chm     =  1.0e9_rp
     comax_chm     = -1.0e9_rp
     kfl_goit2_chm =  0_ip
     alpha1        =  1.0_rp
     alpha2        =  0.0_rp
     alpha3        =  0.0_rp

     !-----------------------------------------------------------------
     !
     ! Allocate memory if necessary
     !
     !-----------------------------------------------------------------

     call chm_explicit_allocate()

     !-----------------------------------------------------------------
     !
     ! Assemble equations:
     !     Gauss-Seidel: Assemble A^i, b^i and solve x^i
     !
     !-----------------------------------------------------------------

     dt(3) = dt(2)
     dt(2) = dt(1)
     dt(1) = 1.0_rp/dtinv

     time_order = min(iter,kfl_tiacc_chm)
     if(      time_order == 1 ) then
        !
        ! 1st order
        !
        alpha1 =  1.0_rp
        alpha2 =  0.0_rp
        alpha3 =  0.0_rp

     else if( time_order == 2 ) then
        !
        ! 2nd order
        !
        dt01   = 1.0_rp
        alpha1 =  ( 1.0_rp + 0.5_rp * dt01)
        alpha2 = -0.5_rp * dt01
        alpha3 =  0.0_rp

     else if( time_order == 3 ) then
        !
        ! 3rd order
        !
        dt012  = dt(1) + dt(2) + dt(3)
        dt01   = dt(1) + dt(2)
        dt12   = dt(2) + dt(3)
        alpha1 = ( 0.25_rp*(dt012*dt01*(dt012+dt01)-dt12*dt(2)*(dt12+dt(2)))&
             &     + 1.0_rp/12.0_rp*(-dt012**3-dt01**3+dt12**3+dt(2)**3))/( dt(2)*dt12)  &
             &   / dt(1)
        alpha2 = ( 0.25_rp*(dt012*dt(1) *(dt012+dt(1))) &
             &     + 1.0_rp/12.0_rp*(-dt012**3-dt(1)**3 +dt12**3))/(-dt(2)*dt(2)) &
             &   / dt(1)
        alpha3 = ( 0.25_rp*(dt01 *dt(1) *(dt01 +dt(1)))&
             &     + 1.0_rp/12.0_rp*(-dt01**3 -dt(1)**3 +dt(2)**3))/(dt(2)*dt12) &
             &   / dt(1)

     end if

     !----------------------------------------------------------------
     !
     ! Matrix assmebly strategy: Monolithic
     !
     !----------------------------------------------------------------

     iclai_chm     = 1
     iclaf_chm     = nclas_chm

     if( INOTMASTER) then

        !
        ! Define unkno(ipoin) from conce(:,iclas,1) for all iclas
        !
        call chm_updunk(ITASK_BEGINN)

        !
        ! Gas phase only
        !
        if (kfl_spray_chm == 0_ip) then
           dt_rho_chm = 0.0_rp
           call chm_matrix()

           call solver_lumped_mass_system(1_ip,dt_rho_chm)

           !
           ! Solve explicit system
           !
           do ipoin = 1,npoin
              do iclas_chm = iclai_chm, iclaf_chm
                 Tim(iclas_chm,ipoin,1) = dt_rho_chm(ipoin)/dt(1)
              end do
           end do

           kpoin = 0
           do ipoin = 1,npoin
              do iclas_chm = iclai_chm, iclaf_chm
                 kpoin             = kpoin + 1
                 rt(iclas_chm,ipoin,1) = bt(kpoin)
              end do
           end do
        !
        ! Spray: liquid + gas phase
        !
        else
           dt_rho_chm = 0.0_rp
           dt_chm     = 0.0_rp
           call chm_matrix()

           call solver_lumped_mass_system(1_ip,dt_rho_chm)
           call solver_lumped_mass_system(1_ip,dt_chm)

           !
           ! Solve explicit system
           !
           do ipoin = 1,npoin
              do iclas_chm = iclai_chm, iclaf_chm-2
                 Tim(iclas_chm,ipoin,1) = dt_rho_chm(ipoin)/dt(1)
              end do
              do iclas_chm = iclaf_chm-1, iclaf_chm
                 Tim(iclas_chm,ipoin,1) = dt_chm(ipoin)/dt(1)
              end do
           end do

           kpoin = 0
           do ipoin = 1,npoin
              do iclas_chm = iclai_chm, iclaf_chm
                 kpoin             = kpoin + 1
                 rt(iclas_chm,ipoin,1) = bt(kpoin)
              end do
           end do
        endif

     end if

     call rhsmod(nclas_chm,rt)

     if( INOTMASTER) then

        reslo = 0.0_rp

        kpoin = 0
        do ipoin = 1,npoin
           do iclas_chm = iclai_chm, iclaf_chm
              kpoin     = kpoin + 1
              reslo     = rt(iclas_chm,ipoin,1) * alpha1 + rt(iclas_chm,ipoin,2) * alpha2 + rt(iclas_chm,ipoin,3) * alpha3
              rholoc    = 0.5_rp*(3.0_rp*Tim(iclas_chm,ipoin,1)-Tim(iclas_chm,ipoin,2))
              tt(kpoin) = tt(kpoin) + ( reslo*rholoc ) * dt(1) / Mass(ipoin)
           end do
        end do
        !
        ! Save old values
        !
        rt (1:nclas_chm,1:npoin,3) = rt (1:nclas_chm,1:npoin,2)
        rt (1:nclas_chm,1:npoin,2) = rt (1:nclas_chm,1:npoin,1)
        Tim(1:nclas_chm,1:npoin,3) = Tim(1:nclas_chm,1:npoin,2)
        Tim(1:nclas_chm,1:npoin,2) = Tim(1:nclas_chm,1:npoin,1)

        !
        ! Dirichet bbcc
        !
        kpoin = 0
        do ipoin = 1,npoin
           do iclas_chm = iclai_chm, iclaf_chm
              kpoin = kpoin + 1
              if( kfl_fixno_chm(iclas_chm,ipoin) > 0 ) &
                   tt(kpoin) = bvess_chm(iclas_chm,ipoin)
           end do
        end do
     end if

     call chm_endite(ITASK_ENDINN)

   end subroutine chm_explicit_solution

 end module mod_chm_explicit

