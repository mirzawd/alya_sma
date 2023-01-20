!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup ChemicMatrixAssembly
!> @{
!> @file    chm_elmope_all.f90
!> @author  Guillaume Houzeaux
!> @brief   Scalar transport element assembly and other element
!>          calculations
!> @details Elemental operations, to include OpenMP, OmpSs,
!>          vectorization and CUDA
!>
!> @}
!------------------------------------------------------------------------
subroutine chm_elmope_all(order)

#include "def_vector_size.inc"
  use def_kintyp,                      only : ip
  use def_domain,                      only : ompss_domain   ! Required by OMPSS
  use def_domain,                      only : lnnod
  use def_domain,                      only : lgaus
#if defined (_OPENMP) || defined (ALYA_OMPSS)
  use mod_parall,                      only : par_omp_nelem_chunk
#endif
  use mod_parall,                      only : num_subd_par
  use mod_parall,                      only : num_pack_par
  use mod_parall,                      only : list_elements_norace_par
  use mod_parall,                      only : num_subd_norace_par
  use mod_parall,                      only : num_pack_norace_par
  use mod_parall,                      only : list_elements_par
  use mod_parall,                      only : typ_list_elements_par
  use def_kermod,                      only : kfl_chemic_vect
  use mod_chm_element_operations,      only : chm_element_operations
  use mod_chm_element_operations_fast, only : chm_element_operations_fast
  use mod_chm_finiteRate,              only : chm_element_operations_finiteRate
  use mod_chm_finiteRate,              only : chm_calc_div_enthalpy_transport_finiteRate
  use mod_chm_finiteRate,              only : chm_calc_hk_grad_Yk_others
  use def_chemic,                      only : kfl_model_chm, kfl_solve_cond_CMC_chm
  use mod_chm_thermophoretic,          only : chm_thermophoretic_calc_nodal_terms
  use def_chemic,                      only : imixf_rk
  use mod_chm_operations_CMC,          only : chm_element_operations_CMC
  use def_chemic,                      only : mixedEq_groups_chm
  use def_chemic,                      only : ngrou_chm

  implicit none

  integer(ip), intent(in)              :: order                     !< What to do
  integer(ip)                          :: isubd,ipack,ielem
  integer(ip)                          :: pnode,pgaus
  integer(ip)                          :: num_subd
  integer(ip)                          :: VECTOR_DIM
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)
  integer(ip)                          :: igrou

  external                             :: runend

  if( order == 1 .or. order == 6 ) then
     !
     ! Element assembly (ORDER=1) and Laplacian assembly (ORDER=6)
     !
     num_subd      =  num_subd_par
     num_pack      => num_pack_par
     list_elements => list_elements_par

  else if( order == 4 ) then
     !
     ! Subgrid scale update
     !
     num_subd      =  num_subd_norace_par
     num_pack      => num_pack_norace_par
     list_elements => list_elements_norace_par

  else

     call runend('CHM_ELMOPE_ALL: NOT CODED')

  end if

  !
  ! Compute gradient of species Yk at the nodes
  !
  if (kfl_model_chm == 3_ip ) &
      call chm_calc_hk_grad_Yk_others

  !
  ! Compute thermophoretic terms at the nodes
  !
  do igrou = 1,ngrou_chm
     if (mixedEq_groups_chm(igrou) % kfl_therm_phor /= 0) &
        call chm_thermophoretic_calc_nodal_terms(mixedEq_groups_chm(igrou) % nequa, &
                                                 mixedEq_groups_chm(igrou) % i_start)
  enddo

  if( order == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Subgrid scale: no race condition
     !
     !-------------------------------------------------------------------

     do isubd = 1,num_subd

        !--------------------------------------------------------------------------
        !$OMP PARALLEL DO                                                         &
        !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                      &
        !$OMP SHARED       ( isubd, par_omp_nelem_chunk )                         &
        !--------------------------------------------------------------------------
        !--------------------------------------------------------------------------
        !$OMP DEFAULT      ( NONE )                                               &
        !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem )                            &
        !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus,kfl_model_chm,kfl_solve_cond_CMC_chm,kfl_chemic_vect )   &
        !--------------------------------------------------------------------------
        !$OMP SHARED       ( order, imixf_rk, VECTOR_DIM                           )
        !--------------------------------------------------------------------------

        do ipack = 1,num_pack(isubd)   !pack = vector

           ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
           pnode = lnnod(ielem)                                ! Number of nodes
           pgaus = lgaus(ielem)                                ! Number of Gauss points

           if( kfl_model_chm == 1  .or. kfl_model_chm == 2_ip ) then
              !
              ! Flamelet model
              !
              if(kfl_chemic_vect == 1_ip) then
                 VECTOR_DIM = size(list_elements(isubd) % packs(ipack) % l,kind=ip)
                 call chm_element_operations_fast(&
                               order,VECTOR_DIM,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)

              else
                 call chm_element_operations(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)
              end if
           else if ( kfl_model_chm == 3 ) then
              !
              ! Finite rate kinetics
              !
              call chm_element_operations_finiteRate(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)

           else if ( kfl_model_chm == 4 ) then
              !
              ! CMC model
              !
              if (kfl_solve_cond_CMC_chm == 1) then
                 call chm_element_operations_CMC(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l, imixf_rk)
              else
                 call chm_element_operations(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)
              end if

           end if

        end do

        !$OMP END PARALLEL DO
     end do

  else

     !-------------------------------------------------------------------
     !
     ! Element assembly and projections: race condition
     !
     !-------------------------------------------------------------------

     do isubd = 1,num_subd

#ifdef ALYA_OMPSS
        num_neigh = size(ompss_domains(isubd) % neighbours)

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
        !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus,kfl_model_chm,kfl_solve_cond_CMC_chm )      &
        !-----------------------------------------------------------------------------
        !$OMP SHARED       ( order, imixf_rk , kfl_chemic_vect, VECTOR_DIM            )
        !-----------------------------------------------------------------------------

        do ipack = 1,num_pack(isubd)

           ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
           pnode = lnnod(ielem)                                ! Number of nodes
           pgaus = lgaus(ielem)                                ! Number of Gauss points

           if( kfl_model_chm == 1  .or. kfl_model_chm == 2_ip ) then
              !
              ! Flamelet model
              !
              if(kfl_chemic_vect == 1_ip) then
                  VECTOR_DIM = size(list_elements(isubd) % packs(ipack) % l,kind=ip)
                  call chm_element_operations_fast(&
                               order,VECTOR_DIM,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)
              else

                  call chm_element_operations(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)
            end if

           else if ( kfl_model_chm == 3 ) then
              !
              ! Finite rate kinetics
              !
              if(kfl_chemic_vect == 1_ip) then

!@                  VECTOR_DIM = size(list_elements(isubd) % packs(ipack) % l,kind=ip)
!@                  call chm_element_operations_finiteRate_fast(&
!@                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)
                   call chm_element_operations_finiteRate(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)


              else
                  call chm_element_operations_finiteRate(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)
              end if
           else if ( kfl_model_chm == 4 ) then
              !
              ! CMC model
              !
              if ( kfl_solve_cond_CMC_chm == 1 ) then
                 call chm_element_operations_CMC(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l, imixf_rk)
              else
                 call chm_element_operations(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)
              end if

           end if

        end do

#ifdef ALYA_OMPSS
        !$OMP END TASK
#else
        !$OMP END PARALLEL DO
#endif
     end do
#ifdef ALYA_OMPSS
     !$OMP  TASKWAIT
#endif
  end if

  !
  ! Compute divergence of enthalpy transport
  !
  if (kfl_model_chm == 3_ip ) &
      call chm_calc_div_enthalpy_transport_finiteRate


end subroutine chm_elmope_all
