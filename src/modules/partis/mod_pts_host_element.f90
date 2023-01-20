!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_host_element.f90
!> @author  houzeaux
!> @date    2020-07-06
!> @brief   Check host element
!> @details Define initial conditions for particle ILAGR
!>          ILAGR: accumulated particle number
!>          KLAGR: local particle number
!>          LAGRTYP(KLAGR) % ILAGR ....... Absolute particle number
!>          LAGRTYP(KLAGR) % ITYPE ....... Particle type
!>          LAGRTYP(KLAGR) % IELEM ....... Host element of particle
!>          LAGRTYP(KLAGR) % KFL_EXIST ... -5 (particle just injected, will be put to 1 right after)
!>          LAGRTYP(KLAGR) % COORD(:) .... Coordinates
!>          LAGRTYP(KLAGR) % DT .......... Current time step dt^n+1
!>          LAGRTYP(KLAGR) % DTO ......... Last time step dt^n
!>          LAGRTYP(KLAGR) % V_FLUID_K ....... Fluid velocity
!>
!>          LAGRTYP(KLAGR) % VELOC ....... Particle velocity at t^n+1
!>
!>          Example with two types:
!>
!>          1. First injection:
!> 
!>          NLAGR_POS =  9
!>          NLAGR_TOT = 18
!>          NLAGR_NEW = 14
!>          o o o   o o o
!>          o o x   o o x
!>          o o x   o o x        
!>          type1   type2
!>
!>          => NLACC_PTS=14
!>
!>          1. Second injection:
!> 
!>          NLAGR_POS =  9
!>          NLAGR_TOT = 18
!>          NLAGR_NEW = 16
!>          x o o   x o o
!>          o o o   o o o
!>          o o o   o o o        
!>          type1   type2
!>
!>          => NLACC_PTS=30
!>
!-----------------------------------------------------------------------

module mod_pts_host_element

  use def_kintyp_basic,            only : ip,rp
  use def_master,                  only : IPARALL
  use def_master,                  only : INOTMASTER
  use def_domain,                  only : ndime
  use def_domain,                  only : mnode
  use def_domain,                  only : meshe
  use def_kermod,                  only : ielse,relse
  use def_kermod,                  only : ndivi
  use mod_memory,                  only : memory_alloca
  use mod_memory,                  only : memory_deallo
  use mod_elsest,                  only : elsest_host_element
  use mod_communications,          only : PAR_POINT_TO_POINT_ARRAY_OPERATION
  use mod_communications,          only : PAR_MAX
  use mod_communications,          only : PAR_SUM
  use def_partis 
  implicit none
  private

  public :: pts_host_element
  public :: pts_host_element_parall
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-06
  !> @brief   Find host elements
  !> @details Find the host elements for a list of particles
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_host_element_parall(&
       nlagr_new,nlagr_inj,particle_place)

    integer(ip),            intent(out)   :: nlagr_new         !< Number of new particles owned by myself
    integer(ip),            intent(inout) :: nlagr_inj         !< New particles
    integer(ip),  pointer,  intent(inout) :: particle_place(:) !< Place of new particles
    integer(ip)                           :: ii,jj
        
    if( IPARALL ) then
       call PAR_POINT_TO_POINT_ARRAY_OPERATION(particle_place, 'IN MY ZONE', 'MIN RANK OR NEGATIVE')
       if (INOTMASTER) then
          nlagr_new = 0
          do ii = 1,nlagr_inj
             
             if(      particle_place(ii) < 0 ) then
                !
                ! One of my neighbors is in charge of this particle
                !
                jj = -particle_place(ii)
                lagrtyp(jj) % kfl_exist = 0
                
             else if( particle_place(ii) > 0 ) then
                !
                ! I take this particle
                !
                nlagr_new = nlagr_new + 1
                
             else if( particle_place(ii) == 0 ) then
                !
                ! Particle not found
                !
                continue
                
             end if
          end do
       end if
       call PAR_MAX(nlagr_inj,'IN MY CODE')
       call PAR_SUM(nlagr_new,'IN MY CODE')
    end if

  end subroutine pts_host_element_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-06
  !> @brief   Find host elements
  !> @details Find the host elements for a list of particles
  !> 
  !-----------------------------------------------------------------------
  
  subroutine pts_host_element(nlagr_pos,nlagr_new,particle_position,host_shapf,host_element)

    integer(ip),           intent(in)    :: nlagr_pos               !< Number of injected particles
    integer(ip),           intent(out)   :: nlagr_new               !< Number of new particles owned by myself
    real(rp),     pointer, intent(in)    :: particle_position(:,:)  !< (ndime,nlagr_pos)
    real(rp),     pointer, intent(inout) :: host_shapf(:,:)
    integer(ip),  pointer, intent(inout) :: host_element(:)
 
    integer(ip)                          :: ielem,ilagr
    real(rp)                             :: coloc(3),deriv(ndime,mnode),xx(3)
    real(rp)                             :: dummr(3),relse_sav,dista

    !----------------------------------------------------------------------
    ! 
    ! Look for host elements
    !
    !----------------------------------------------------------------------

    relse_sav   = relse(1)
    !relse(1)    = 0.0_rp
    nlagr_new   = 0
    xx          = 0.0_rp

    !$OMP PARALLEL DO SCHEDULE (DYNAMIC,500)                          & 
    !$OMP DEFAULT   (NONE)                                            &    
    !$OMP PRIVATE   (ilagr,xx,dummr,ielem,deriv,coloc,dista)          &
    !$OMP SHARED    (nlagr_pos,particle_position,ielse,relse,         &
    !$OMP            kfl_exacs_pts,meshe,ndivi,                       &
#ifndef NDIMEPAR
    !$OMP            ndime,                                           &
#endif
    !$OMP            host_shapf,host_element)                         &
    !$OMP REDUCTION (+:nlagr_new)
    do ilagr = 1,nlagr_pos
       xx(1:ndime) = particle_position(1:ndime,ilagr)
       if( kfl_exacs_pts /= 0 ) then
          call pts_exacso(1_ip,0.0_rp,dummr,dummr,xx)
       end if
       call elsest_host_element(&
            ielse,relse,1_ip,meshe(ndivi),xx,ielem,&
            host_shapf(:,ilagr),deriv,coloc,dista)
       if( ielem > 0 ) then
          host_element(ilagr) = ielem
          nlagr_new           = nlagr_new + 1 ! Number of new particles injected successfully
       end if
    end do
    !$OMP END PARALLEL DO

    relse(1)  = relse_sav

  end subroutine pts_host_element

end module mod_pts_host_element
!> @}
