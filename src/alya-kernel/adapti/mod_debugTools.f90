!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_debugTools.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   Simple module to allow/remove outputing in adaptivity
!> @details Simple module to allow/remove outputing in adaptivity
!-----------------------------------------------------------------------
MODULE mod_debugTools
  !
  use def_kintyp_basic,       only: ip,rp,lg
  !
  implicit none
  
  logical(lg), parameter :: out_debug_text                = .false.!.true.!
  logical(lg), parameter :: out_debug_paraview            = .false.!.true.!
  logical(lg), parameter :: out_debug_paraview_firstLast  = .false.!.false.!.true.!
  integer(ip)            :: out_id_firstLast = 0_ip
  !
  ! Performance variables:
  !
  logical(lg), parameter :: out_performance      = .false.
  logical(lg), parameter :: do_out_large_profile = .false.
  !
  real(rp) :: time_total, time_split, time_coll, time_swap, time_smooth ! time per iter
  real(rp) :: deb_timeTotal        = 0_rp
  real(rp) :: deb_collTimeTotal    = 0_rp
  real(rp) :: deb_splitTimeTotal   = 0_rp
  real(rp) :: deb_swapTimeTotal    = 0_rp
  real(rp) :: deb_optiTimeTotal    = 0_rp
  !
  real(rp) :: deb_setNewMeshModify = 0.0_rp
  real(rp) :: deb_metricMerge      = 0.0_rp
  real(rp) :: deb_qualityMesh      = 0.0_rp
  real(rp) :: deb_minQCavity       = 0.0_rp
  real(rp) :: deb_minQSubset       = 0.0_rp
  real(rp) :: deb_qualSubset       = 0.0_rp
  real(rp) :: deb_distSubset       = 0.0_rp
  real(rp) :: deb_nodeToElem       = 0.0_rp
  real(rp) :: deb_isBouFromIni     = 0.0_rp
  real(rp) :: deb_optiFun          = 0.0_rp
  real(rp) :: deb_BLS              = 0.0_rp
  real(rp) :: deb_FunMin           = 0.0_rp
  real(rp) :: deb_numer            = 0.0_rp
  real(rp) :: deb_swapLoop         = 0.0_rp
  real(rp) :: deb_swapRest         = 0.0_rp
  real(rp) :: deb_swapIni          = 0.0_rp
  real(rp) :: deb_tagEdgesByQ      = 0.0_rp
  real(rp) :: deb_swapPreTagQ      = 0.0_rp
  real(rp) :: deb_swapPreTagNE     = 0.0_rp
  real(rp) :: deb_swapSetMesh      = 0.0_rp
  real(rp) :: deb_splitLoop        = 0.0_rp
  real(rp) :: deb_splitTagEdges    = 0.0_rp
  real(rp) :: deb_splitEdgeLen     = 0.0_rp
  real(rp) :: deb_collLoop         = 0.0_rp
  real(rp) :: deb_collTagEdges     = 0.0_rp
  real(rp) :: deb_collEdgeLen      = 0.0_rp
  
  real(rp) :: deb_meshEdgeLen      = 0.0_rp
  real(rp) :: deb_meshEdges        = 0.0_rp
  real(rp) :: deb_allEdgesLen      = 0.0_rp
  !
  !
  !
  public
  !
  !
  !
CONTAINS
  !
  !
  !
  subroutine out_performance_profile()
    implicit none
  
    if(out_performance) then
     write(*,1) time_total,(time_coll+time_split+time_swap+time_smooth),time_coll,time_split, time_swap, time_smooth
     1  format(/,&
     5x,'|---------------------------|',/, &
     5x,'| Total time     : ',f8.2,' |',/, &
     5x,'|   Time_oper    : ',f8.2,' |',/, &
     5x,'|     time_coll  : ',f8.2,' |',/, &
     5x,'|     time_split : ',f8.2,' |',/, &
     5x,'|     time_swap  : ',f8.2,' |',/, &
     5x,'|     time_smooth: ',f8.2,' |',/, &
     5x,'|---------------------------|',/&
     )
    end if
  
  end subroutine out_performance_profile
  !
  !
  !
  subroutine out_performance_profile_global()
    implicit none
    
  
    if(out_performance) then
     write(*,1) deb_timeTotal,&
                deb_collTimeTotal,  deb_splitTimeTotal, deb_swapTimeTotal, deb_optiTimeTotal
     1  format(/,&
     5x,'|---------------------------|',/, &
     5x,'| Total        : ',f8.2,'   |',/, &
     5x,'|---------------------------|',/, &
     5x,'| Global time per operation |',/, &
     5x,'|   time_coll  : ',f8.2,'   |',/, &
     5x,'|   time_split : ',f8.2,'   |',/, &
     5x,'|   time_swap  : ',f8.2,'   |',/, &
     5x,'|   time_smooth: ',f8.2,'   |',/, &
     5x,'|---------------------------|')
   
     if(do_out_large_profile) then
       
       write(*,2) deb_setNewMeshModify, deb_metricMerge, &   ! mesh and metric modification operations
                  deb_nodeToElem,deb_isBouFromIni,&          ! topological calculation (neighbors, boundaries)
                  deb_minQSubset, deb_qualSubset, deb_qualityMesh, deb_minQCavity , &  ! quality computation
                  deb_optiFun, deb_distSubset, deb_numer, deb_FunMin, deb_BLS, &!   mesh optimization (opti and distortion)
                  deb_swapLoop, deb_splitLoop, deb_collLoop
       2  format(&
       5x,'| setNewMeshMod: ',f8.2,'   |',/, &
       5x,'| metric%merge : ',f8.2,'   |',/, &
       5x,'| nodeToElem   : ',f8.2,'   |',/, &
       5x,'| isBouFromIni : ',f8.2,'   |',/, &
       5x,'| minQSubset   : ',f8.2,'   |',/, &
       5x,'| qualSubset   : ',f8.2,'   |',/, &
       5x,'| mesh quality : ',f8.2,'   | -> rellevant',/, &
       5x,'| minQCavity   : ',f8.2,'   | -> rellevant',/, &
       5x,'| optiFunction : ',f8.2,'   | -> rellevant (smoothing)',/, &
       5x,'|   distSubset : ',f8.2,'   | -> rellevant (smoothing)',/, &
       5x,'|   DerNumer   : ',f8.2,'   |',/, &
       5x,'|   FMin(NR/SD): ',f8.2,'   |',/, &
       5x,'|   BackLS     : ',f8.2,'   |',/, &
       5x,'| swapLoop     : ',f8.2,'   |',/, &
       5x,'| splitLoop    : ',f8.2,'   |',/, &
       5x,'| collLoop     : ',f8.2,'   |',/, &
       5x,'|---------------------------|')
     
       write(*,3) deb_collTimeTotal,  deb_collLoop,   deb_collTagEdges,  deb_collEdgeLen, &
                  deb_splitTimeTotal, deb_splitLoop, deb_splitTagEdges, deb_splitEdgeLen, &
                  deb_swapTimeTotal,                 &
                  deb_swapPreTagQ,                   &
                  deb_swapIni,  deb_tagEdgesByQ,     &
                  deb_swapLoop, deb_swapSetMesh,     &
                  deb_swapPreTagNE!,                  &
                  !deb_swapRest
                
       3  format(&
       5x,'| collTimeTotal: ',f8.3,'   |',/, &
       5x,'|   collLoop   : ',f8.3,'   |',/, &
       5x,'|   tagEdges   : ',f8.3,'   |',/, &
       5x,'|     edgeLen  : ',f8.3,'   |',/, &
       5x,'| splitTimeTot : ',f8.3,'   |',/, &
       5x,'|   splitLoop  : ',f8.3,'   |',/, &
       5x,'|   tagEdges   : ',f8.3,'   |',/, &
       5x,'|     edgeLen  : ',f8.3,'   |',/, &
       5x,'| swapTimeTotal: ',f8.3,'   |',/, &
       5x,'|   preTag_Q   : ',f8.3,'   |',/, &
       5x,'|   swapIni    : ',f8.3,'   |',/, &
       5x,'|     tagByQ   : ',f8.3,'   |',/, &
       5x,'|   swapLoop   : ',f8.3,'   |',/, &
       5x,'|   setMesh    : ',f8.3,'   |',/, &
       5x,'|   nodeToElem : ',f8.3,'   |',/, &!     5x,'|   swapRest   : ',f8.2,'   |',/, &
       5x,'|---------------------------|')
     
     
       write(*,4) deb_meshEdgeLen,  deb_meshEdges, deb_allEdgesLen
       4  format(&
       5x,'|---------------------------|',/, &
       5x,'| Mesh Edge Lengths         |',/, &
       5x,'|   Total      : ',f8.3,'   |',/, &
       5x,'|      edges   : ',f8.3,'   |',/, &
       5x,'|      lengths : ',f8.3,'   |',/, &
       5x,'|---------------------------|')
       
     end if
     
    end if
  
  end subroutine out_performance_profile_global
  !
  !
  !
END MODULE mod_debugTools