!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_minle_time_step.f90
!> @author  Gerard Guillamet
!> @brief   Critical time step based on min. lenght (Vectorized version)
!> @details Critical time step based on min. lenght (Vectorized version)
!> @}
!-----------------------------------------------------------------------

module mod_sld_minle_time_step

#include "def_vector_size.inc"
  use def_kintyp,              only : ip,rp
  use def_domain,              only : ndime,ntens,mnode,elmar,npoin,lnods
  use def_domain,              only : vmass,mnode,coord,lgaus,ltype,nnode
  use def_domain,              only : nelem,lnnod,lelch
  use def_domain,              only : ompss_domain   ! Required by OMPSS
  use def_domain,              only : ompss_domains
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  use mod_element_integration, only : element_shape_function_derivatives_jacobian_vector_2
  use mod_elmgeo_vector,       only : elmgeo_cartesian_derivatives_jacobian_vectorized_cpu
  
  use mod_parall,              only : par_omp_nelem_chunk
  use mod_parall,              only : num_subd_par
  use mod_parall,              only : num_pack_par
  use mod_parall,              only : list_elements_norace_par
  use mod_parall,              only : num_subd_norace_par
  use mod_parall,              only : num_pack_norace_par
  use mod_parall,              only : list_elements_par
  use mod_parall,              only : typ_list_elements_par
  
  implicit none

  private

  public :: sld_minle_time_step_all

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-07-27
  !> @brief   Global time step subroutine based on minimum lenght 
  !> @details Global time step subroutine based on minimum lenght
  !>
  !-----------------------------------------------------------------------

  subroutine sld_minle_time_step_all(dtmin)

    use def_elmtyp,                  only : ELFEM, ELINT
    use mod_communications,          only : PAR_MIN
    use def_solidz,                  only : celen_sld
    use mod_sld_interface_element,   only : ELINT_stable_time_increment
    
    real(rp),    intent(inout)           :: dtmin
    integer(ip)                          :: isubd,ipack,ielem,kelem
    integer(ip)                          :: pnode
    integer(ip)                          :: VECTOR_DIM
    integer(ip)                          :: num_subd
#ifdef ALYA_OMPSS
    integer(ip)                          :: num_neigh
#endif
    integer(ip),                 pointer :: num_pack(:)
    type(typ_list_elements_par), pointer :: list_elements(:)

    real(rp)                             :: dtcri
    
    num_subd      =  num_subd_par
    num_pack      => num_pack_par
    list_elements => list_elements_par

    dtmin = huge(1.0_rp)
    
    do isubd = 1,num_subd

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
       !$OMP PRIVATE      ( ipack,pnode,ielem,VECTOR_DIM )                          &
       !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus )
       !-----------------------------------------------------------------------------

       do ipack = 1,num_pack(isubd)

          ielem      = list_elements(isubd) % packs(ipack) % l(1)                 ! Select first element
          pnode      = lnnod(ielem)                                               ! Number of nodes
          VECTOR_DIM = size(list_elements(isubd) % packs(ipack) % l,kind=ip)

          if(      lelch(ielem) == ELFEM ) then
#if defined PNODE_VALUE
             call sld_minle_time_step_element_operations_ELFEM(&
                  VECTOR_DIM,&
                  list_elements(isubd) % packs(ipack) % l,dtcri)
#else
             call sld_minle_time_step_element_operations_ELFEM(&
                  VECTOR_DIM,pnode,&
                  list_elements(isubd) % packs(ipack) % l,dtcri)        
#endif
          else if( lelch(ielem) == ELINT ) then
             do kelem = 1, size(list_elements(isubd) % packs(ipack) % l)
                ielem = list_elements(isubd) % packs(ipack) % l(kelem)
                if( ielem > 0_ip ) then
                   call ELINT_stable_time_increment(ielem,dtcri)
                   celen_sld(ielem) = 1.0_rp ! TODO: for interface element should take the min. in-plane size.
                end if
             end do
          end if
          !
          ! Min. time increment
          !
          dtmin = min( dtmin, dtcri )

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
    ! Parall: Look for minimum over all subdomains (dtmin)
    !
    call PAR_MIN(dtmin,'IN MY CODE')
 
  end subroutine sld_minle_time_step_all
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-09-02
  !> @brief   Element operations for minimum element length strategy
  !> @details Element operations for minimum element length strategy
  !> 
  !-----------------------------------------------------------------------

#if defined PNODE_VALUE
  subroutine sld_minle_time_step_element_operations_ELFEM(VECTOR_DIM,list_elements,dtcri)
#else
  subroutine sld_minle_time_step_element_operations_ELFEM(VECTOR_DIM,pnode,list_elements,dtcri)  
#endif    
    use def_master,                  only : ITER_K
    use def_master,                  only : displ
    use def_domain,                  only : lnods, ltype, lmate, lorde
    use def_domain,                  only : coord, hnatu
    use mod_communications,          only : PAR_MIN
    use mod_elmgeo,                  only : elmgeo_element_characteristic_length
    use def_solidz,                  only : celen_sld,velas_sld,kfl_celen_sld
    use mod_sld_vect_assembly_elfem, only : sld_vect_element_length
    
    implicit none
    
    integer(ip),          intent(in)  :: VECTOR_DIM
#if defined PNODE_VALUE                        
    integer(ip), parameter            :: pnode = PNODE_VALUE           !< Number of nodes
#else
    integer(ip),          intent(in)  :: pnode                         !< Number of nodes
#endif
    integer(ip), pointer, intent(in)  :: list_elements(:)
    real(rp),             intent(out) :: dtcri

    integer(ip)                       :: ielem0(1)
    integer(ip)                       :: list_elements_p(VECTOR_DIM)   ! List of elements (always positive)
    !
    ! Gather
    !
    real(rp)                          :: elcod(VECTOR_DIM,ndime,mnode)
    !
    ! Indices and dimensions
    !
    integer(ip)                       :: ielem,ivect,inode,ipoin
    integer(ip)                       :: pelty,porde,pmate
    real(rp)                          :: hleng(VECTOR_DIM,ndime)
    real(rp)                          :: dt

#define DEF_VECT 1:VECTOR_SIZE_CPU
    !
    ! Initialization
    !
    ielem = list_elements(1)
    pelty = ltype(ielem)     
    porde = lorde(pelty)     
    pmate = lmate(ielem)   
    !
    ! List of elements. Put last non-zero element in the list to
    ! avoid zero arrays
    !
    list_elements_p = list_elements
    ielem0 = minloc(list_elements,list_elements==0)
    if( ielem0(1) > 1 ) list_elements_p(ielem0(1):VECTOR_DIM) = list_elements(ielem0(1)-1)

    !--------------------------------------------------------------------
    !
    ! Gather
    !
    !--------------------------------------------------------------------

    if( kfl_celen_sld == 1 ) then
       do ivect = 1,VECTOR_DIM
          ielem = list_elements_p(ivect)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,ITER_K)
          end do
       end do
    else
       do ivect = 1,VECTOR_DIM
          ielem = list_elements_p(ivect)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)
          end do
       end do
    end if

    !--------------------------------------------------------------------
    !
    ! Element length
    !
    !--------------------------------------------------------------------

    call elmgeo_element_characteristic_length(&
         ndime,pnode,elmar(pelty)%dercg(:,:),elcod(:,:,:),hleng(:,:),hnatu(pelty))

    !--------------------------------------------------------------------
    !
    ! Min. length and critical time step
    !
    !--------------------------------------------------------------------

    dtcri = huge(1.0_rp) 
    do ivect = 1,VECTOR_DIM
       ielem = list_elements_p(ivect)
       celen_sld(ielem) = hleng(ivect,ndime) / real(porde,rp)
       dt = celen_sld(ielem)/velas_sld(1,pmate)
       dtcri = min(dtcri,dt)
    end do
    
  end subroutine sld_minle_time_step_element_operations_ELFEM

end module mod_sld_minle_time_step
