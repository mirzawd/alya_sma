!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_matrix.f90
!> @author  Solidz Team
!> @date    November, 2017-Adds doxygen
!> @brief   This routine computes the matrix and right hand side
!> @details
!>      OUTPUT
!>      USES
!>        sld_elmope
!>        sld_bouope
!>
!>      USED BY
!>        mod_solution_methods_sld
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_matrix(itask)

   use def_kintyp,             only : ip, rp, lg
   use def_master,             only : INOTMASTER, coupling, IMASTER,ITER_K
   use def_master,             only : npoi1, npoi2, npoi3, mem_modul
   use def_master,             only : rhsid, amatr, INOTEMPTY
   use def_master,             only : gdedet, gdeinv
   use def_master,             only : ITASK_MATRIX, solve_sol, modul
   use def_master,             only : CPU_ASSEMBLY,CPU_MINI_ASSEMBLY
   use def_master,             only : CPU_MAXI_ASSEMBLY,CPU_COUNT_ASSEMBLY
   use def_master,             only : cutim, timei, timef
   use def_master,             only : gdepo
   use def_domain,             only : ndime, npoin
   use def_domain,             only : vmass
   use mod_memory,             only : memory_alloca
   use mod_memory,             only : memory_deallo, memory_alloca_min
   use def_solver,             only : SOL_MATRIX_HAS_CHANGED
   use def_solidz,             only : kfl_penal_sld, kfl_limit_sld, kfl_vofor_sld, kfl_conta_sld
   use def_solidz,             only : fexte_sld, finte_sld, macce_sld, fcont_sld, SLD_PDN_UNILATERAL
   use def_solidz,             only : fext2_sld, fint2_sld, fine2_sld
   use def_solidz,             only : kfl_gdepo
   use def_solidz,             only : vofor_sld, cpu_ass_sol_sld
   use def_solidz,             only : ndofn_sld, frxid_sld
   use def_solidz,             only : nzrhs_sld, vmasx_sld, vmass_sld
   use def_solidz,             only : fnatu_sld
   use def_solidz,             only : kfl_vecto_sld
   use def_solidz,             only : kfl_timet_sld
   use def_solidz,             only : SLD_EXPLICIT_SCHEME, SLD_IMPLICIT_SCHEME  
   use def_solidz,             only : pressEndoIntegralOld_sld, pressEndoIntegralNew_sld
   use def_solidz,             only : pushForwardIntegralOld_sld, pushForwardIntegralNew_sld
   use mod_ker_timeline,       only : ker_timeline
   use mod_communications,     only : PAR_BARRIER, PAR_SUM
   use mod_timings,            only : timings_assembly
   use mod_alya2dlb,           only : alya2dlb_DLB_Enable
   use mod_alya2dlb,           only : alya2dlb_DLB_Disable
   use mod_array_operations,   only : array_operations_norm2
   use mod_array_operations,   only : array_operations_initialization
   use mod_eccoupling,         only : exm_sld_ecc_initialize_troponin, &
   &                                  exm_sld_ecc_interchange_troponin, &
   &                                  kfl_exmsld_ecc

   implicit none

   integer(ip), intent(in)         :: itask
   integer(ip)                     :: izrhs,ipoin,idime,itott,ierr
   real(rp)                        :: time1,time2,time3,time4
   real(rp),    pointer            :: frnat(:,:)

   call ker_timeline('INI_ASSEMBLY')

   !-----------------------------------------------------------------
   !
   ! Initializations
   !
   !-----------------------------------------------------------------

   time1     =  0.0_rp
   time2     =  0.0_rp
   time3     =  0.0_rp

   if( INOTMASTER ) then
      !
      ! Initialize solver
      !
      if (itask == 1_ip) then        ! explicit

         do izrhs=1,nzrhs_sld
            rhsid(izrhs)=0.0_rp
         end do
         call inisol()

      else if(itask == 2_ip) then    ! implicit

         do izrhs=1,nzrhs_sld
            rhsid(izrhs)=0.0_rp
         end do
         call inisol()

      end if

   end if

   !
   ! Initialize arrays
   !
   if( itask <= 2_ip .and. INOTEMPTY ) then
      !
      ! Mass matrices
      !
      call array_operations_initialization(vmass_sld) ! Diagonal mass matrix
      call array_operations_initialization(vmasx_sld) ! Lumped mass matrix (X-FEM)
      !
      ! Force vectors
      !
      call array_operations_initialization(frxid_sld) ! Reaction forces
      call array_operations_initialization(finte_sld) ! Internal force vector
      call array_operations_initialization(fexte_sld) ! External force vector
      call array_operations_initialization(macce_sld) ! Inertial force vector
      call array_operations_initialization(fcont_sld) ! Contact force (reaction) vector
      !
      ! Global gradient deformation operator
      !
      if( kfl_gdepo /= 0 ) then
         call array_operations_initialization(gdepo)
         call array_operations_initialization(gdeinv)
         call array_operations_initialization(gdedet)
      end if

   end if

   !-----------------------------------------------------------------
   !
   ! EXMEDI coupling : initialize variables to be assembled
   !
   !-----------------------------------------------------------------

   if( kfl_exmsld_ecc )then
      if( itask <= 2 .and. INOTEMPTY )then
         call exm_sld_ecc_initialize_troponin()
      endif
   end if

   !-----------------------------------------------------------------
   !
   ! Matrix, preconditioner, RHS and elment assemblies
   !
   !-----------------------------------------------------------------

   if( INOTMASTER ) ierr = alya2dlb_DLB_Enable()

   call cputim(time1)
   
   if( itask <= 2 ) then
      if( kfl_vecto_sld ) then
         if(      kfl_timet_sld == SLD_EXPLICIT_SCHEME ) then
            call sld_vect_elmope_expl(itask)
         else if( kfl_timet_sld == SLD_IMPLICIT_SCHEME ) then
            call sld_vect_elmope_impl(itask)
         end if
      else
         call sld_elmope(itask)
      end if
   end if
   
   call cputim(time2)

#ifdef ALYA_DLB
   call PAR_BARRIER()
   if( INOTMASTER ) ierr = alya2dlb_DLB_Disable()
#endif

   !-----------------------------------------------------------------
   !
   ! Solve push forward (slaves)
   !
   !-----------------------------------------------------------------

   if( INOTEMPTY ) then

      if( kfl_gdepo /= 0 ) then
         call rhsmod(ndime*ndime,gdepo)
         call rhsmod(ndime*ndime,gdeinv)
         call rhsmod(       1_ip,gdedet)
         do ipoin=1, npoin
            gdepo(1:ndime,1:ndime,ipoin)  = gdepo( 1:ndime,1:ndime,ipoin) / vmass(ipoin)
            gdeinv(1:ndime,1:ndime,ipoin) = gdeinv(1:ndime,1:ndime,ipoin) / vmass(ipoin)
            gdedet(ipoin)                 = gdedet(ipoin) / vmass(ipoin)
         end do
      end if

   end if

   !-----------------------------------------------------------------
   !
   ! EXMEDI coupling : interchange assembled variables
   !
   !-----------------------------------------------------------------

   if( coupling('SOLIDZ','EXMEDI') >= 1_ip .or. coupling('EXMEDI','SOLIDZ') >= 1_ip ) then
      if( itask <= 2 .and. INOTEMPTY )then
         call exm_sld_ecc_interchange_troponin()
      endif
   end if

   !-----------------------------------------------------------------
   !
   ! Boundary assembly
   !
   !-----------------------------------------------------------------

   call cputim(time3)
   if( INOTMASTER ) then
      call sld_bouope(itask)
   end if
   ! Sum slaves' contributions and update pressure surface integrals
   call PAR_SUM(ndime,pressEndoIntegralNew_sld(1:ndime),'IN MY CODE')
   call PAR_SUM(      pushForwardIntegralNew_sld       ,'IN MY CODE')
   pressEndoIntegralOld_sld(1:3)   = pressEndoIntegralNew_sld(1:3) 
   pushForwardIntegralOld_sld      = pushForwardIntegralNew_sld

   call cputim(time4)

   !-----------------------------------------------------------------
   !
   ! Other stuff
   !
   !-----------------------------------------------------------------

   if( INOTMASTER ) then
      !
      ! Coupling
      !
      call sld_coupli(ITASK_MATRIX)
      !
      ! Limiter
      !
      if (kfl_limit_sld==1) call sld_limite()
      !
      ! Add external volume forces
      !
      if( kfl_vofor_sld > 0 ) then
         !
         ! Internal nodes
         !
         do ipoin = 1,npoin
            if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
               itott = (ipoin-1) * ndofn_sld
               do idime = 1,ndime
                  itott = itott + 1
                  rhsid(itott) = rhsid(itott) + vofor_sld(idime,ipoin)*(cutim/(timef-timei))
                  fexte_sld(itott) = fexte_sld(itott) + vofor_sld(idime,ipoin)*(cutim/(timef-timei))
               end do
            end if
         end do

      end if
      !
      ! Penalize system
      !
      if (kfl_penal_sld == 1) call sld_pensys(amatr)

   end if
   !
   ! Parallel exchage: rhsmod for force vectors and mass matrix
   !
   if( INOTEMPTY ) then
      call rhsmod(ndime,frxid_sld)
      call rhsmod(ndime,finte_sld)
      call rhsmod(ndime,fexte_sld)
      call rhsmod(ndime,macce_sld)
      call rhsmod(1_ip, vmass_sld)
      if ( kfl_conta_sld == SLD_PDN_UNILATERAL ) call rhsmod(ndime, fcont_sld)
   end if
   !
   ! Sum
   !
   !
   ! L2 norms of the force vectors required for the converge criteria
   !
   fint2_sld = array_operations_norm2(finte_sld,DOFS=ndime)
   fext2_sld = array_operations_norm2(fexte_sld,DOFS=ndime)
   fine2_sld = array_operations_norm2(macce_sld,DOFS=ndime)
   !
   ! Timings
   !
   cpu_ass_sol_sld(1) = time2 - time1
   cpu_ass_sol_sld(4) = time4 - time3
   call timings_assembly(cpu_ass_sol_sld(1),cpu_ass_sol_sld(4))
   !
   ! L2 norm of natural force imposed in solver
   !
   nullify(frnat)
   if( solve_sol(1) % kfl_bvnat == 1_ip .and. associated(solve_sol(1) % bvnat) ) then
      if( IMASTER ) then
         call memory_alloca_min(mem_modul(1:2,modul),'FRNAT','sld_matrix',frnat)
      else
         call memory_alloca(mem_modul(1:2,modul),'FRNAT','sld_matrix',frnat,ndime,npoin,'DO_NOT_INITIALIZE')
      end if
      frnat = 0.0_rp
      fnatu_sld = 0.0_rp
      if( INOTMASTER ) then
         frnat = solve_sol(1) % bvnat
         call rhsmod(ndime,frnat)
      end if
      fnatu_sld = array_operations_norm2(frnat,DOFS=ndime)
      call memory_deallo(mem_modul(1:2,modul),'FRNAT','sld_matrix',frnat)
   else
      fnatu_sld = 0.0_rp
   end if
   !
   ! Timeline
   !
   call ker_timeline('END_ASSEMBLY')
   !
   ! Matrix has been assembled
   !
   solve_sol(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED ! Matrix has been assembled

end subroutine sld_matrix
