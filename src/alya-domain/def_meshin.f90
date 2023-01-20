!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_meshin

  ! NAME
  !   def_domain
  ! DESCRIPTION
  !   This module is the header of the meshing part
  !***
  !-----------------------------------------------------------------------
  use def_kintyp

  !------------------------------------------------------------------------
  ! Connectivity arrays
  !------------------------------------------------------------------------
  integer(ip), pointer              ::     ptoel1(:)=>null()
  integer(ip), pointer              ::     ptoel2(:)=>null()
  integer(ip), pointer              ::     ptosi1(:)=>null(),ptosi2(:)=>null()
  integer(ip), pointer              ::     ptosi1old(:)=>null(),ptosi2old(:)=>null()
  integer(ip), pointer              ::     eltoel(:,:)=>null(),ledge(:,:)=>null()
  integer(ip), pointer              ::     sitosi(:,:)=>null(),lftoed(:,:)=>null()
  integer(ip), pointer              ::     lstof(:)=>null(),lsold(:,:)=>null() 
  !------------------------------------------------------------------------
  ! Boundary layer mesh arrays
  !------------------------------------------------------------------------
  integer(ip), pointer              ::     ldiag(:)=>null()
  real(rp), pointer                 ::     rblay(:)=>null() 
  !------------------------------------------------------------------------
  ! Geometry 
  !------------------------------------------------------------------------
  real(rp), pointer                 ::     rnopo(:,:)=>null(),rnofa(:,:)=>null()
  integer(ip),pointer               ::     lboup(:,:)=>null(),lsurf(:)=>null()
  integer(ip),pointer               ::     lpsur(:)=>null(),lpsid(:)=>null()
  integer(ip),pointer               ::     lline(:)=>null()
  integer(ip),pointer               ::     lside(:,:),lptype(:,:)  
  integer(ip)                       ::     nsurf,nboup 
  integer(8)                        ::     memor_msh(2)  

  !------------------------------------------------------------------------
  ! Cartesian mesh
  !------------------------------------------------------------------------
  type(cell), pointer               ::     lcell(:)=>null()
  integer(ip)                       ::     ncell,kfl_ifbox,mleve,nsour
  integer(ip)                       ::     nvoxx,nvoxy,nvoxz
  integer(ip),pointer               ::     btoel1(:)=>null(),btoel2(:)=>null()
  real(rp)                          ::     bboxbin(3,2)
  real(rp)                          ::     rsuni,rscal,boxin(3,2)
  !------------------------------------------------------------------------
  ! Surface mesh
  !------------------------------------------------------------------------
  integer(ip),pointer               ::     lface(:,:)=>null() 
  integer(ip),pointer               ::     lfloc(:,:)=>null(),lsurfloc(:)=>null() 
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------
  integer(ip),parameter             ::     ID_CORNER=4         
  integer(ip),parameter             ::     ID_RIDGE=3         
  integer(ip),parameter             ::     ID_CUSP=2         
  integer(ip),parameter             ::     ID_SMOOTH=1         
  !------------------------------------------------------------------------
  ! Volume mesh
  !------------------------------------------------------------------------
  integer(ip),pointer               ::     lcart(:)=>null()
  integer(ip),pointer               ::     lptri(:)=>null()
  integer(ip),pointer               ::     lelem(:)=>null()
  integer(ip),pointer               ::     lmark(:)=>null()
  integer(ip),pointer               ::     elem(:,:)=>null()
  real(rp),pointer                  ::     rsize(:)=>null()
  real(rp),pointer                  ::     coor(:,:)=>null()
  !------------------------------------------------------------------------
  ! Frontal approach
  !------------------------------------------------------------------------
  integer(ip),pointer               ::     lfront(:,:)=>null()
  real(rp),pointer                  ::     rfront(:)=>null()
  integer(ip),pointer               ::     lheap(:)=>null()
  integer(ip),pointer               ::     lfmark(:)=>null()
  integer(ip),pointer               ::     lfapo(:,:)=>null()
  integer(ip),pointer               ::     lpofa(:)=>null()
  integer(ip),pointer               ::     lfhole(:)=>null()
  integer(ip),pointer               ::     lfahol(:)=>null()
  !------------------------------------------------------------------------
  ! Sources
  !------------------------------------------------------------------------
  real(rp),pointer                  ::     rsgeo(:,:,:)=>null()
  real(rp),pointer                  ::     rsour(:,:)=>null()




end module def_meshin
