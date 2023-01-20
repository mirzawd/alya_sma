!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_updtss.f90
!> @author  Guillaume Houzeaux
!> @brief   Critical time steo
!> @details Computes the time step size for the incompressible NS
!>          equation.
!> @}
!-----------------------------------------------------------------------
subroutine nsi_updtss()
  use def_kintyp,              only : ip,rp
  use def_master,              only : INOTMASTER
  use def_master,              only : kfl_timco,dtinv,veloc
  use def_master,              only : veloc_forw
  use def_master,              only : times
  use def_kermod,              only : kfl_adj_prob
  use def_domain,              only : lnods,coord,nelem
  use def_domain,              only : ltype,nnode,elmar,ngaus
  use def_domain,              only : mnode,lmate,ndime
  use def_domain,              only : lorde,hnatu
  use def_nastin,              only : kfl_timei_nsi
  use def_nastin,              only : kfl_advec_nsi
  use def_nastin,              only : kfl_ellen_nsi
  use def_nastin,              only : dtcri_nsi
  use def_nastin,              only : dtsgs_nsi,kfl_dttyp_nsi
  use def_nastin,              only : kfl_sgsti_nsi,dtinv_nsi
  use def_nastin,              only : safet_nsi,tamin_nsi
#if defined(_OPENMP)
  use mod_parall,              only : par_omp_nelem_chunk
#endif
  use mod_nsi_eigen_time_step, only : nsi_eigen_time_step_all
  use mod_communications,      only : PAR_MIN,PAR_BARRIER
  use mod_alya2dlb,            only : alya2dlb_DLB_Enable
  use mod_alya2dlb,            only : alya2dlb_DLB_Disable
  implicit none
  
  integer(ip)      :: ielem,inode,ipoin,ierr
  integer(ip)      :: pnode,pelty,pmate,pgaus
  real(rp)         :: dtcri
  real(rp)         :: dtmin
  real(rp)         :: cartd(ndime,mnode)
  real(rp)         :: chave(ndime,2)
  real(rp)         :: elcod(ndime,mnode)
  real(rp)         :: elvel(ndime,mnode)
  real(rp)         :: chale(2)
  real(rp)         :: hleng(3)
  real(rp)         :: tragl(9)
  integer(ip)      :: porde
  
  if( INOTMASTER ) ierr = alya2dlb_DLB_Enable()
#ifdef AVOID_DTC_CALC
  dtcri_nsi = 1.0
  return
#endif

  call times(7) % ini()

  if( kfl_timei_nsi /= 0 ) then

     if( kfl_dttyp_nsi == 1 ) then

        !---------------------------------------------------------------
        !
        ! Minimum of element time steps
        !
        !---------------------------------------------------------------

        if( INOTMASTER ) then
           
           dtmin = huge(1.0_rp)

           !--------------------------------------------------------------------------
           !$OMP PARALLEL  DO                                                        &
           !$OMP SCHEDULE  ( DYNAMIC , par_omp_nelem_chunk )                         &
           !$OMP DEFAULT   ( NONE )                                                  &
           !$OMP PRIVATE   ( ielem,pelty,pnode,porde,pmate,inode,ipoin,pgaus,        &
           !$OMP             elcod,elvel,tragl,hleng,chave,chale,cartd,dtcri )       &
           !$OMP SHARED    ( ltype,nnode,lorde,kfl_ellen_nsi,kfl_adj_prob,           &
           !$OMP             lmate,lnods,coord,veloc,elmar,hnatu,kfl_advec_nsi,      &
           !$OMP             veloc_forw,nelem,ngaus,                                 &
#ifndef NDIMEPAR
           !$OMP             ndime,                                                  &
#endif
           !$OMP             par_omp_nelem_chunk )                                   &
           !$OMP REDUCTION ( MIN:dtmin )
           !--------------------------------------------------------------------------
           !
           ! Loop over elements
           !
           elements: do ielem = 1,nelem
              pelty = ltype(ielem)

              if( pelty > 0 ) then
                 pnode = nnode(pelty)
                 porde = lorde(pelty)
                 pgaus = ngaus(pelty)
                 pmate = lmate(ielem)
                 !
                 ! Default: ELVEL and ELCOD
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                    if (kfl_adj_prob == 0) then
                       elvel(1:ndime,inode) = veloc(1:ndime,ipoin,1)
                    else
                       elvel(1:ndime,inode) = veloc_forw(1:ndime,ipoin,1)
                    endif
                 end do
                 !
                 ! HLENG and TRAGL at center of gravity
                 !
                 call elmlen(&
                      ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                      hnatu(pelty),hleng)
                 !
                 ! Compute the characteristic length: CHALE
                 !
                 call elmchl(&
                      tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                      porde,hnatu(pelty),kfl_advec_nsi,kfl_ellen_nsi)
                 !
                 ! Time step
                 !
                 call nsi_elmtss(&
                      pelty,pmate,pnode,lnods(1,ielem),ielem,elcod,elvel,&
                      cartd,chale,hleng,dtcri)

                 dtmin = min( dtmin , dtcri )

              end if
           end do elements
           !$OMP END PARALLEL DO

        end if
        !
        ! Parall: Look for minimum over all subdomains (dtmin)
        !
        call PAR_MIN(dtmin,'IN MY CODE')
        
     else if( kfl_dttyp_nsi == 2 ) then

        !---------------------------------------------------------------
        !
        ! Adaptive time step
        !
        !---------------------------------------------------------------

        call nsi_updtse(dtmin)

     else if( kfl_dttyp_nsi == 3 ) then

        !---------------------------------------------------------------
        !
        ! Use minimum TAU
        !
        !---------------------------------------------------------------

        dtmin = tamin_nsi

     else if( kfl_dttyp_nsi == 4 ) then

        !---------------------------------------------------------------
        !
        ! Eigenvalues
        !
        !---------------------------------------------------------------

        call nsi_eigen_time_step_all(dtmin)
        
     end if
     !
     ! Assign 1/dt
     !
     dtcri_nsi = dtmin
     if( dtcri_nsi /= 0.0_rp ) dtinv_nsi = 1.0_rp/(dtcri_nsi*safet_nsi)
     dtsgs_nsi = dtinv_nsi
     if( kfl_sgsti_nsi == 0 ) dtsgs_nsi = 0.0_rp
     if( kfl_timco     == 1 ) dtinv = max(dtinv,dtinv_nsi)

  end if

#ifdef ALYA_DLB
  call PAR_BARRIER()
    if( INOTMASTER ) ierr = alya2dlb_DLB_Disable()
#endif
    
  call times(7) % add()

    call times(7) % add()

end subroutine nsi_updtss
 
