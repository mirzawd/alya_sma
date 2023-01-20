!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_energy.f90
!> @author  Gerard Guillamet
!> @date    February, 2018
!>           - Subroutine written
!> @brief   Calculation of energies and energetic balance checks
!> @details Calculation of energies and energetic balance checks
!>
!>          \verbatim
!>          The following energies are implemented:
!>           - ALLIE: All internal energy
!>           - ALLWK: All external work
!>           - ALLKE: All kinetic energy
!>           - ETOTA: Total energy
!>
!>          Energy balance:
!>           - Energetic balance according to Belytschko
!>          \endverbatim
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures.
!>
!> @todo    To do list::\n
!>           - Energetic balance for subdomains (locally, see Belytschko)
!>           - Include energy from natural forces (fsi problems)
!>           - Include strain energy
!>           - Include artificial energy from rayleigh damping
!>
!> @}
!------------------------------------------------------------------------

module mod_sld_energy

  use def_kintyp,         only : ip, rp
  use def_domain,         only : ndime, npoin
  use def_master,         only : INOTEMPTY, ITER_K, TIME_N
  use def_solidz,         only : fintt_sld, fextt_sld
  use def_solidz,         only : allie_sld, allwk_sld, allke_sld, etota_sld

  implicit none

  private

  public :: sld_energy
  public :: sld_updene

contains

  !-----------------------------------------------------------------------
  !>
  !> @brief   This routine calculates energies and perform a an energy
  !>          balance check.
  !> @details
  !>          \verbatim
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine sld_energy()

    use def_master,         only : ittim,rhsid
    use mod_operations,     only : operations_parallel_dot_product
    use def_master,         only : modul, mem_modul
    use def_master,         only : times
    use mod_memory,         only : memory_alloca
    use mod_memory,         only : memory_deallo
    use def_solidz,         only : SLD_DYNAMIC_PROBLEM
    use def_solidz,         only : kfl_timei_sld
    use def_solidz,         only : veloc_sld, ddisp_sld
    use def_solidz,         only : vmass_sld
    use def_solidz,         only : eener_sld
    use def_solidz,         only : ndofn_sld

    implicit none

    integer(ip)                 :: ipoin,idime,idofn
    real(rp)                    :: allwk, allie, allke          ! Energies
    real(rp)                    :: eener
    real(rp), pointer           :: fintsum(:,:)
    real(rp), pointer           :: fextsum(:,:)
    real(rp), pointer           :: rhstota(:,:)
    real(rp), pointer           :: velosqu(:,:)
    real(rp), pointer           :: vmassto(:,:)
    real(rp), pointer           :: ddisp_aux(:,:)
    
    call times(2) % ini()

    nullify(fintsum)
    nullify(fextsum)
    nullify(rhstota)
    nullify(velosqu)
    nullify(vmassto)
    nullify(ddisp_aux)
    !
    ! Allocate memory
    !
    call memory_alloca(mem_modul(1:2,modul),'FINTSUM','sld_energy',fintsum,ndime,npoin)
    call memory_alloca(mem_modul(1:2,modul),'FEXTSUM','sld_energy',fextsum,ndime,npoin)
    call memory_alloca(mem_modul(1:2,modul),'RHSTOTA','sld_energy',rhstota,ndime,npoin)
    call memory_alloca(mem_modul(1:2,modul),'VELOSQU','sld_energy',velosqu,ndime,npoin)
    call memory_alloca(mem_modul(1:2,modul),'VMASSTO','sld_energy',vmassto,ndime,npoin)
    ddisp_aux => ddisp_sld(:,:,ITER_K) ! To be consistent with the interface of the parallel dot product
    !
    ! Preliminary calculations
    !
    if( INOTEMPTY ) then
       !
       ! Forces summation between current and previous iteration and rhs in array form
       do ipoin = 1,npoin
          idofn = (ipoin-1)*ndofn_sld
          fintsum(1:ndime,ipoin) = fintt_sld(1:ndime,ipoin,ITER_K) + fintt_sld(1:ndime,ipoin,TIME_N)
          fextsum(1:ndime,ipoin) = fextt_sld(1:ndime,ipoin,ITER_K) + fextt_sld(1:ndime,ipoin,TIME_N)
          ! rhs in array form
          rhstota(1:ndime,ipoin) = rhsid(idofn+1:idofn+ndime)
       end do
       ! Square of the velocity and lumped mass matrix
       if (kfl_timei_sld == SLD_DYNAMIC_PROBLEM) then
          ! Square of velocity and lumped mass matrix
          do ipoin = 1,npoin
             do idime = 1,ndime
                velosqu(idime,ipoin) = veloc_sld(idime,ipoin,ITER_K)*veloc_sld(idime,ipoin,ITER_K)
                vmassto(idime,ipoin) = vmass_sld(ipoin)
             end do
          end do
       end if

    end if

    !
    ! All energies (whole model)
    !
    ! All internal energy
    call operations_parallel_dot_product(fintsum,ddisp_aux,allie,'IN MY CODE')

    if (ittim == 1_ip) then
       allie_sld(ITER_K) = allie_sld(TIME_N) + abs(allie)
    else
       allie_sld(ITER_K) = allie_sld(TIME_N) + 0.5_rp*abs(allie)
    end if
    ! All work (external energy)
    call operations_parallel_dot_product(fextsum,ddisp_aux,allwk,'IN MY CODE')
    if (ittim == 1_ip) then
       allwk_sld(ITER_K) = allwk_sld(TIME_N) + abs(allwk)
    else
       allwk_sld(ITER_K) = allwk_sld(TIME_N) + 0.5_rp*abs(allwk)
    end if
    ! All kinetic energy
    if (kfl_timei_sld == SLD_DYNAMIC_PROBLEM) then
       call operations_parallel_dot_product(velosqu,vmassto,allke,'IN MY CODE')
       allke_sld(ITER_K) = 0.5_rp*allke
    end if

    !
    ! Total energy
    !
    etota_sld(ITER_K) = allie_sld(ITER_K) + allke_sld(ITER_K) - allwk_sld(ITER_K)

    !
    ! Error in energy flow
    !
    call operations_parallel_dot_product(rhstota,ddisp_aux,eener,'IN MY CODE')
    eener_sld = abs(eener)

    !
    ! Energy balance
    !
    if (etota_sld(ITER_K) <= 1.0e-2_rp*max(allie_sld(ITER_K),allwk_sld(ITER_K),allke_sld(ITER_K)) ) then
       !print*,'Energy balance OK'
    end if
    !
    ! Deallocate memory of local arrays
    !
    call memory_deallo(mem_modul(1:2,modul),'FINTSUM','sld_energy',fintsum)
    call memory_deallo(mem_modul(1:2,modul),'FEXTSUM','sld_energy',fextsum)
    call memory_deallo(mem_modul(1:2,modul),'RHSTOTA','sld_energy',rhstota)
    call memory_deallo(mem_modul(1:2,modul),'VELOSQU','sld_energy',velosqu)
    call memory_deallo(mem_modul(1:2,modul),'VMASSTO','sld_energy',vmassto)

    call times(2) % add()

  end subroutine sld_energy

  !-----------------------------------------------------------------------
  !>
  !> @brief   This routine performs several updates for energies and
  !>          variables related to them.
  !> @details
  !>          \verbatim
  !>          ITASK = ITASK_BEGSTE ... Initial guess outer iterations
  !>                  ITASK_BEGITE ... Initial guess inner iterations
  !>                  ITASK_ENDITE ... Updates (end interation)
  !>                  ITASK_ENDINN ... Updates (end iteration converged)
  !>                  ITASK_ENDSTE ... End time step
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine sld_updene(itask)

    use def_master,           only : ITER_AUX
    use def_master,           only : ITASK_BEGSTE, ITASK_ENDSTE
    use def_master,           only : ITASK_BEGITE, ITASK_ENDITE
    use def_solidz,           only : nprev_sld

    implicit none

    integer(ip), intent(in)       :: itask !< Update variables at selected case

    select case (itask)

    case( ITASK_BEGSTE )

       !------------------------------------------------------------------
       !
       !  Outer iterations (begin time step)
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTEMPTY) then
          fintt_sld(1:ndime,1:npoin,ITER_AUX) = fintt_sld(1:ndime,1:npoin,nprev_sld)
          fextt_sld(1:ndime,1:npoin,ITER_AUX) = fextt_sld(1:ndime,1:npoin,nprev_sld)
       end if
       ! Energies
       allie_sld(ITER_AUX) = allie_sld(nprev_sld)
       allwk_sld(ITER_AUX) = allwk_sld(nprev_sld)
       allke_sld(ITER_AUX) = allke_sld(nprev_sld)
       etota_sld(ITER_AUX) = etota_sld(nprev_sld)

    case( ITASK_BEGITE )

       !------------------------------------------------------------------
       !
       !  Inner iterations (begin iterations)
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTEMPTY) then
          fintt_sld(1:ndime,1:npoin,ITER_K) = fintt_sld(1:ndime,1:npoin,ITER_AUX)
          fextt_sld(1:ndime,1:npoin,ITER_K) = fextt_sld(1:ndime,1:npoin,ITER_AUX)
       end if
       ! Energies
       allie_sld(ITER_K) = allie_sld(ITER_AUX)
       allwk_sld(ITER_K) = allwk_sld(ITER_AUX)
       allke_sld(ITER_K) = allke_sld(ITER_AUX)
       etota_sld(ITER_K) = etota_sld(ITER_AUX)

    case( ITASK_ENDITE )

       !------------------------------------------------------------------
       !
       !  Updates (end iterations converged)
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTEMPTY) then
          fintt_sld(1:ndime,1:npoin,ITER_AUX) = fintt_sld(1:ndime,1:npoin,ITER_K)
          fextt_sld(1:ndime,1:npoin,ITER_AUX) = fextt_sld(1:ndime,1:npoin,ITER_K)
       end if
       ! Energies
       allie_sld(ITER_AUX) = allie_sld(ITER_K)
       allwk_sld(ITER_AUX) = allwk_sld(ITER_K)
       allke_sld(ITER_AUX) = allke_sld(ITER_K)
       etota_sld(ITER_AUX) = etota_sld(ITER_K)

    case( ITASK_ENDSTE )

       !------------------------------------------------------------------
       !
       !  End Time Step
       !
       !------------------------------------------------------------------
       ! Force vectors
       if (INOTEMPTY) then
          fintt_sld(1:ndime,1:npoin,TIME_N) = fintt_sld(1:ndime,1:npoin,ITER_K)
          fextt_sld(1:ndime,1:npoin,TIME_N) = fextt_sld(1:ndime,1:npoin,ITER_K)
       end if
       ! Energies
       allie_sld(TIME_N) = allie_sld(ITER_K)
       allwk_sld(TIME_N) = allwk_sld(ITER_K)
       allke_sld(TIME_N) = allke_sld(ITER_K)
       etota_sld(TIME_N) = etota_sld(ITER_K)

    end select

  end subroutine sld_updene

end module mod_sld_energy

