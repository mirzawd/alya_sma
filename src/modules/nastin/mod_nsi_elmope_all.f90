!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_nsi_elmope_all 

  implicit none
  public :: nsi_elmope_all 

contains
!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_elmope_all.f90
!> @author  Guillaume Houzeaux
!> @brief   Navier-Stokes system element assembly and other element
!>          calculations
!> @details Elemental operations, to include OpenMP, OmpSs,
!>          vectorization and CUDA
!>
!> @}
!------------------------------------------------------------------------
  subroutine nsi_elmope_all(itask)
    use def_kintyp,                 only : ip,rp
    use def_master,                 only : IEMPTY
    use def_domain,                 only : ompss_domain   ! Required by OMPSS
#if !defined(OPENACCHHH) && defined(ALYA_OMPSS)
    use def_domain,                 only : ompss_domains
#endif
    use def_domain,                 only : lnnod
    use def_domain,                 only : lgaus
#ifdef OPENACCHHH
    use def_domain,                 only : ndime
#endif
    use def_nastin,                 only : resis_nsi
    use def_nastin,                 only : itsta_nsi
    use def_nastin,                 only : resgs_nsi
    use def_nastin,                 only : tamin_nsi
    use def_nastin,                 only : rmsgs_nsi
    use def_nastin,                 only : tamax_nsi
    use def_nastin,                 only : dtmax_nsi
    use def_nastin,                 only : cputi_assembly_nsi
    use def_nastin,                 only : kfl_assem_nsi
#if defined(_OPENMP)
    use mod_parall,                 only : par_omp_nelem_chunk
#endif
    use mod_parall,                 only : num_subd_par
    use mod_parall,                 only : num_pack_par
    use mod_parall,                 only : list_elements_norace_par
    use mod_parall,                 only : num_subd_norace_par
    use mod_parall,                 only : num_pack_norace_par
    use mod_parall,                 only : list_elements_par
    use mod_parall,                 only : typ_list_elements_par

    use mod_communications,                  only : PAR_BARRIER
    use mod_nsi_element_operations,          only : nsi_element_operations
    use mod_nsi_element_operations_fast,     only : nsi_element_operations_fast
    use mod_nsi_element_operations_test,     only : nsi_element_operations_fast8
    use def_master,                          only : velom
    use def_parall,                          only : kfl_openacc_streams_per_dev
    use def_domain,                          only : npoin
    use def_kermod,                          only : kfl_noslw_ker
#ifdef _OPENACC
    use def_master,                          only : veloc
    use def_master,                          only : rhsid
    use def_kermod,                          only : kfl_nswel_ker
    use def_nastin,                          only : dt_rho_nsi
    use def_nastin,                          only : mass_rho_nsi
    use def_nastin,                          only : ncomp_nsi
#endif

#ifndef NASTIN_PRIVATE_OFF
    use mod_nsi_element_operations_fast_dev, only : nsi_element_operations_fast_dev
    use mod_nsi_element_operations_hh71,     only : nsi_element_operations_hh71
    use mod_nsi_element_operations_hh80,     only : nsi_element_operations_hh80
    use mod_nsi_element_operations_hh90,     only : nsi_element_operations_hh90
    use mod_nsi_element_operations_hh91,     only : nsi_element_operations_hh91
#endif


    implicit none

    integer(ip), intent(in)              :: itask                     !< What to do
    integer(ip)                          :: isubd,ipack,ielem
    integer(ip)                          :: pnode,pgaus
#if !defined(OPENACCHHH) && defined(ALYA_OMPSS)
    integer(ip)                          :: jsubd
    integer(ip)                          :: num_neigh
#endif
    integer(ip)                          :: num_subd
    integer(ip),                 pointer :: num_pack(:)
    type(typ_list_elements_par), pointer :: list_elements(:)
    real(rp)                             :: time_detail(10)
#ifdef OPENACCHHH
    integer(ip)                          :: streamid
#endif
    integer(ip)                          :: npoin_min

#if defined USE_LIKWID
    call likwid_markerStartRegion("elmope")
#endif

    if( IEMPTY ) return

    dtmax_nsi =-1.0_rp

    npoin_min = max(1_ip, npoin)

    if( itask == 1 .or. itask == 6 ) then
       !
       ! Element assembly (ITASK=1) and Laplacian assembly (ITASK=6)
       !
       num_subd      =  num_subd_par
       num_pack      => num_pack_par
       list_elements => list_elements_par

    else if( itask == 4 ) then
       !
       ! Subgrid scale update
       !
       num_subd      =  num_subd_norace_par
       num_pack      => num_pack_norace_par
       list_elements => list_elements_norace_par

    else

       call runend('NSI_ELMOPE_ALL: NOT CODED')

    end if

    time_detail = 0.0_rp

    if( itask == 4 ) then

       !-------------------------------------------------------------------
       !
       ! Subgrid scale: no race condition
       !
       !-------------------------------------------------------------------

       do isubd = 1,num_subd    !colors

          !--------------------------------------------------------------------------
          !$OMP PARALLEL DO                                                         &
          !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                      &
          !$OMP SHARED       ( isubd, par_omp_nelem_chunk )                         &
          !--------------------------------------------------------------------------
          !--------------------------------------------------------------------------
          !$OMP DEFAULT      ( NONE )                                               &
          !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem )                            &
          !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus )                 &
          !--------------------------------------------------------------------------
          !$OMP SHARED       ( itask )                                              &
          !$OMP REDUCTION    ( +:resis_nsi,itsta_nsi,resgs_nsi,time_detail )        &
          !$OMP REDUCTION    ( MIN:tamin_nsi )                                      &
          !$OMP REDUCTION    ( MAX:rmsgs_nsi,tamax_nsi,dtmax_nsi                    )
          !--------------------------------------------------------------------------

          do ipack = 1,num_pack(isubd)   !pack = vector

             ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
             pnode = lnnod(ielem)                                ! Number of nodes
             pgaus = lgaus(ielem)                                ! Number of Gauss points

             call nsi_element_operations(&
                  itask,pnode,pgaus,list_elements(isubd) % packs(ipack) % l,&
                  resis_nsi,itsta_nsi,resgs_nsi,tamin_nsi,rmsgs_nsi,tamax_nsi,&
                  dtmax_nsi,time_detail)

          end do

          !$OMP END PARALLEL DO
       end do

    else !itask /= 4

       !-------------------------------------------------------------------
       !
       ! Element assembly and projections: race condition
       !
       !-------------------------------------------------------------------

       !$acc enter data copyin(rhsid,veloc(1:ndime,1:npoin,1:1))
       !$acc update device(rhsid,veloc(1:ndime,1:npoin,1:1))

       if ( kfl_assem_nsi > 20 ) then
          !$acc enter data copyin(dt_rho_nsi(1:npoin_min))
          !$acc enter data copyin(mass_rho_nsi(1:npoin,1:ncomp_nsi))
          !$acc update device(dt_rho_nsi,mass_rho_nsi)
       end if

       if(associated(velom)) then
!!! !$acc enter data create(veloc(1:ndime,1:npoin))
!!! !$acc update device(veloc(1:ndime,1:npoin))
       end if



       if ( kfl_noslw_ker /= 0_ip ) then
!!! !!$acc enter data copyin (kfl_nswel_ker(1:npoin))
       end if
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
          !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus,kfl_assem_nsi )      &
          !-----------------------------------------------------------------------------
          !$OMP SHARED       ( itask )                                                 &
          !$OMP REDUCTION    ( +:resis_nsi,itsta_nsi,resgs_nsi,time_detail )           &
          !$OMP REDUCTION    ( MIN:tamin_nsi )                                         &
          !$OMP REDUCTION    ( MAX:rmsgs_nsi,tamax_nsi,dtmax_nsi                       )
          !-----------------------------------------------------------------------------
#endif

          do ipack = 1,num_pack(isubd)

             ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
             pnode = lnnod(ielem)                                ! Number of nodes
             pgaus = lgaus(ielem)                                ! Number of Gauss points

             select case ( kfl_assem_nsi )

             case ( 5_ip ) 
                !
                ! Assembly 5 (dev)
                !
                if(kfl_openacc_streams_per_dev /= 1) then

#ifdef OPENACCHHH

#ifdef _OPENACC
                   streamid=(mod(ipack,kfl_openacc_streams_per_dev) + 1)
#endif
                   call nsi_element_operations_fast8(&
                        size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                        pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail,streamid)
#else
#ifndef NASTIN_PRIVATE_OFF              
                   call nsi_element_operations_fast_dev(&
                        size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                        pnode,pgaus,list_elements(isubd) % packs(ipack) %l,time_detail)
#else
                   call runend("nastin++ sources missing")
#endif
#endif
                else
#ifndef NASTIN_PRIVATE_OFF              
                   call nsi_element_operations_fast_dev(&
                        size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                        pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
#else
                   call runend("nastin++ sources missing")
#endif
                end if

             case ( 2_ip )
                !
                ! Assembly 2 (stable)
                !
                call nsi_element_operations_fast(&
                     size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                     pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
#ifndef NASTIN_PRIVATE_OFF      

             case ( 71_ip )
                !
                ! Assembly 71
                !
                call nsi_element_operations_hh71(&
                     size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                     pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
             
             case ( 80_ip )
                !
                ! Assembly 80
                !
                call nsi_element_operations_hh80(&
                     size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                     pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
             
             case ( 90_ip )
                !
                ! Assembly 90
                !
                call nsi_element_operations_hh90(&
                     size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                     pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)

             case ( 91_ip )
                !
                ! Assembly 91
                !            
                call nsi_element_operations_hh91(&
                     size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                     pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
#else
             case ( 71 )
                !
                ! Assembly 71
                !
                call runend("nastin++ sources missing")
#endif
             case default
                !
                ! Assembly by default
                !
                call nsi_element_operations(&
                     itask,pnode,pgaus,list_elements(isubd) % packs(ipack) % l,&
                     resis_nsi,itsta_nsi,resgs_nsi,tamin_nsi,rmsgs_nsi,tamax_nsi,&
                     dtmax_nsi,time_detail)

             end select

          end do

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
          !$OMP END TASK
#else
          !$OMP END PARALLEL DO
#endif
#endif

       end do
#if defined USE_LIKWID
       call   likwid_markerStopRegion("elmope")
#endif

#ifndef OPENACCHHH


#ifdef ALYA_OMPSS
       !$OMP  TASKWAIT
#endif
#else

       !$acc wait 

#endif
       !
       ! Accumulate CPU time
       !
       cputi_assembly_nsi = cputi_assembly_nsi + time_detail

    end if

    !$acc update host(rhsid)
    if ( kfl_assem_nsi > 20 ) then
       !$acc update host(dt_rho_nsi,mass_rho_nsi)
    end if

!!!  !!$acc exit data delete(dt_rho_nsi(1:npoin_min),mass_rho_nsi(1:npoin,1:ncomp_nsi), & 
!!!!  !!$acc                  rhsid,veloc(1:ndime,1:npoin,1:ncomp_nsi))

    if ( kfl_noslw_ker /= 0_ip ) then
!!! !!$acc exit data delete (kfl_nswel_ker(1:npoin))
    end if

    if( associated(velom)) then
!!! !$acc exit data delete(veloc(1:ndime,1:npoin))
    end if

  end subroutine nsi_elmope_all
 
end module mod_nsi_elmope_all
!> @}
