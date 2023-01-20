!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_elmope(itask)

  use def_kintyp,                 only : ip,rp
  use def_domain,                 only : lnnod
  use def_domain,                 only : lgaus
#ifdef ALYA_OMPSS
  use def_domain,                 only : ompss_domain   ! Required by OMPSS
  use def_domain,                 only : ompss_domains
#elif defined(_OPENMP)
  use mod_parall,                 only : par_omp_nelem_chunk
#endif
  use mod_parall,                 only : num_subd_par
  use mod_parall,                 only : num_pack_par
  use mod_parall,                 only : list_elements_par
  use mod_parall,                 only : typ_list_elements_par
  use mod_communications,         only : PAR_BARRIER

  implicit none

  integer(ip), intent(in)              :: itask                     !< What to do
  integer(ip)                          :: isubd,ipack,ielem,kelem
  integer(ip)                          :: pnode,pgaus

  integer(ip)                          :: ielem_min
  real(rp)                             :: gpdet_min

#ifdef ALYA_OMPSS
  integer(ip)                          :: jsubd
  integer(ip)                          :: num_neigh
#endif
  integer(ip)                          :: num_subd
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)
  real(rp)                             :: time_detail(10)

  num_subd      =  num_subd_par
  num_pack      => num_pack_par
  list_elements => list_elements_par

  time_detail = 0.0_rp
  gpdet_min   = huge(1.0_rp)
  ielem_min   = 0

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
     !$OMP DEFAULT      ( NONE )                                                  &
     !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem,kelem )                         &
     !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus )                    &
     !-----------------------------------------------------------------------------
     !$OMP SHARED       ( itask,ielem_min,gpdet_min )
     !-----------------------------------------------------------------------------

     do ipack = 1,num_pack(isubd)

        ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
        pnode = lnnod(ielem)                                ! Number of nodes
        pgaus = lgaus(ielem)                                ! Number of Gauss points

        do kelem = 1,size(list_elements(isubd) % packs(ipack) % l)
           ielem = list_elements(isubd) % packs(ipack) % l(kelem)
           if( ielem > 0 ) &
                call sld_element_operations(&
                itask,ielem,ielem_min,gpdet_min)
        end do

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
  !
  ! Detection of Negative jacobian
  !
  if( ielem_min > 0 ) call sld_detjac(ielem_min,gpdet_min)

end subroutine sld_elmope
