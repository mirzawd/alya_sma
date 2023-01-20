!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_bouope_all.f90
!> @author  Guillaume Houzeaux
!> @brief   Matrix assembly: boundary contribution
!> @details Boundary operations
!> @}
!-----------------------------------------------------------------------
subroutine nsi_bouope_all(itask)
  
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin 
  use mod_ker_proper
#if !defined(OPENACCHHH) && defined(_OPENMP)
  use mod_parall,                       only : par_omp_nboun_chunk
#endif

  use mod_parall,                       only : num_subd_nboun_par
  use mod_parall,                       only : num_pack_nboun_par
  use mod_parall,                       only : list_boundaries_par
  use mod_parall,                       only : typ_list_elements_par
  use mod_nsi_boundary_operations,      only : nsi_boundary_operations
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast5
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast_bck1
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast_bck2
  use mod_memory
  use mod_nsi_frixgb      ! machine learning
  implicit none 
  integer(ip), intent(in)              :: itask
  integer(ip)                          :: isubd,ipack,iboun
  integer(ip)                          :: pnodb,pgaub
#if !defined(OPENACCHHH) && defined(ALYA_OMPSS)
  integer(ip)                          :: jsubd
  integer(ip)                          :: num_neigh
#endif
  integer(ip)                          :: pblty,ivect
  integer(ip)                          :: num_subd
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_boundaries(:)
  real(rp)                             :: vikin
  integer(ip)                          :: dummi
  real(rp),                allocatable :: tvenos(:),local_re_all(:),dim_u_all(:)
  real(rp),                allocatable :: ustar_all(:)
! this is a bad idea
  real(rp)                :: gbden(4)
  real(rp)                :: gbvis(4)

  !
  ! Loop indices
  !
  num_subd        =  num_subd_nboun_par
  num_pack        => num_pack_nboun_par
  list_boundaries => list_boundaries_par

  if( IMASTER ) return
  !
  ! Initialize nodal traction and boundary mass
  !
  if( itask == 1_ip ) then
     if( associated(notra_nsi) ) notra_nsi = 0.0_rp
     if( associated(massb_nsi) ) massb_nsi = 0.0_rp
  end if
  !
  ! Loop over boundaries
  !
  do isubd = 1,num_subd
  ! machine learning
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (kfl_mlwm_ker == 1_ip) then
       if (.not.allocated(tveno_all)) allocate(tveno_all(mgaub, nboun))
       if (.not.allocated(tvscd_all)) allocate(tvscd_all(mgaub, nboun))
       if (.not.allocated(dim_u)) allocate(dim_u(mgaub, nboun))
       if (.not.allocated(local_re)) allocate(local_re(mgaub, nboun))
       if (.not.allocated(delta_ml)) allocate(delta_ml(nboun*mgaub))
       if (.not.allocated(delsc_ml)) allocate(delsc_ml(nboun*mgaub))
       if (.not.allocated(dim_y)) allocate(dim_y(nboun*mgaub))

!      Initialise:
       tveno_all = 0.0_rp
       tvscd_all = 0.0_rp
       dim_u     = 0.0_rp
       local_re  = 0.0_rp
       delta_ml  = 0.0_rp
       delsc_ml  = 0.0_rp
       dim_y     = 1.0_rp
       do ipack = 1,num_pack(isubd)
          iboun = list_boundaries(isubd) % packs(ipack) % l(1)                     ! Select first element
          pnodb = lnnob(iboun)                                                     ! Number of boundary nodes
          pblty = ltypb(iboun)                                                     ! Type of boundary
          pgaub = ngaus(pblty)                                                     ! Number of Gauss points
          if (iboun > 0_ip)  then
           call ker_proper('DENSI','PGAUB',dummi,iboun,gbden,pnodb,pgaub,elmar(pblty)%shape)
           call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis,pnodb,pgaub,elmar(pblty)%shape)       
           vikin = gbvis(1)/gbden(1)
           call nsi_ml_ustar_all(&
               size(list_boundaries(isubd) % packs(ipack) % l,KIND=ip),&
               list_boundaries(isubd) % packs(ipack) % l,&
               vikin,pnodb,pgaub) 
          end if
       end do
          tvenos = reshape(tveno_all,[size(tveno_all)])
          local_re_all = reshape(local_re,[size(tveno_all)])
          dim_u_all = reshape(dim_u,[size(tveno_all)])
          allocate(ustar_all(size(tvenos)))
          ustar_all = 0.0_rp
!          if (size(delta_ml) .ne. 0_ip) call frixgb(local_re_all,dim_y,dim_u_all,ustar_all,size(dim_y),kfl_paral)  
          if (size(delta_ml) .ne. 0_ip) call frixgb(local_re,dim_y,dim_u,ustar_all,size(dim_y))  
          ustars = reshape(ustar_all,[size(tveno_all,1),size(tveno_all,2)])
          deallocate(local_re_all)
          deallocate(dim_u_all)
          deallocate(delta_ml)
          deallocate(delsc_ml)
          deallocate(tveno_all)
          deallocate(tvscd_all)
          deallocate(dim_y)
     end if

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
     num_neigh = size(ompss_boundaries(isubd) % neighbours,KIND=ip)

     !-----------------------------------------------------------------------------------
     !$OMP TASK         COMMUTATIVE(                                                    &
     !$OMP              [ompss_boundaries(ompss_boundaries(isubd) % neighbours(jsubd)), &
     !$OMP              jsubd = 1,num_neigh] ) PRIORITY(num_neigh)                      &
     !$OMP FIRSTPRIVATE ( num_neigh,jsubd,isubd )                                       &
     !$OMP SHARED       ( ompss_boundaries )                                            &
     !-----------------------------------------------------------------------------------
#else 
     !-----------------------------------------------------------------------------------
     !$OMP PARALLEL DO                                                                  &
     !$OMP SCHEDULE     ( DYNAMIC , par_omp_nboun_chunk )                               &
     !$OMP SHARED       ( isubd,par_omp_nboun_chunk )                                   &
     !-----------------------------------------------------------------------------------
#endif
     !-----------------------------------------------------------------------------------
     !$OMP DEFAULT      ( SHARED )                                                      &
     !$OMP PRIVATE      ( ipack,pnodb,pgaub,iboun,pblty,ivect )                         &
     !$OMP SHARED       ( list_boundaries,num_pack,ltypb,lnnob,ngaus,lelbo)             &
     !-----------------------------------------------------------------------------------
     !$OMP SHARED       ( itask                                                         )
     !-----------------------------------------------------------------------------------
#endif

     do ipack = 1,num_pack(isubd)

        iboun = list_boundaries(isubd) % packs(ipack) % l(1)                     ! Select first element
        pnodb = lnnob(iboun)                                                     ! Number of boundary nodes
        pblty = ltypb(iboun)                                                     ! Type of boundary
        pgaub = ngaus(pblty)                                                     ! Number of Gauss points

#ifndef AVOID_BOUN_OPER        
        if ( kfl_asbou_nsi == 5 ) then
           call nsi_boundary_operations_fast5(&
                itask,size(list_boundaries(isubd) % packs(ipack) % l,KIND=ip),pnodb,&
                pgaub,list_boundaries(isubd) % packs(ipack) % l)           
        else
           do ivect = 1,size(list_boundaries(isubd) % packs(ipack) % l,KIND=ip)
              iboun = list_boundaries(isubd) % packs(ipack) % l(ivect)           ! Boundary
              if( iboun > 0 ) &
                   call nsi_boundary_operations(&
                   itask,pnodb,pgaub,list_boundaries(isubd) % packs(ipack) % l(ivect:))
           end do
        end if
#endif
        
     end do
!     if (allocated(ustar_all)) deallocate(ustar_all)
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
#endif
  
end subroutine nsi_bouope_all
