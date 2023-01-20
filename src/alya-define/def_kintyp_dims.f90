!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup Kinds_and_types
!> Kinds ands types of Alya
!> @{
!> @file    def_kintyp.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   Definition of dimensions
!> @details Definition of mesh dimensions
!>
!-----------------------------------------------------------------------

module def_kintyp_dims

  use def_kintyp_basic,  only : ip
  use def_elmtyp,        only : element_max
  
  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip), parameter  :: nelty = element_max ! # of element types
  integer(ip), parameter  :: mfiel = 500         ! Maximum number of fields

  !------------------------------------------------------------------------
  !
  ! Dimensions: read in readim
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_autbo,             &      ! Automatic boundaries
       npoin,                 &      ! # of nodal points
       nelem,                 &      ! # of elements
       necnt,                 &      ! # of contact elements
       nncnt,                 &      ! # of contact nodes
       nboun,                 &      ! # of boundary elements
       nperi,                 &      ! # periodic nodes
       lexis(nelty),          &      ! List of existing elements
       nexis,                 &      ! Number of different element types in the geometry file
       utype,                 &      ! Value of the unique type, -1 if several types
       npoib,                 &      ! # immersed nodes
       nboib,                 &      ! # immersed boundaries
       lexib(nelty),          &      ! List of existing IB elements
       nhang,                 &      ! # hanging nodes
       nimbo,                 &      ! # IB
       nrbod,                 &      ! # RB
       nzone,                 &      ! Number of zones
       nsubd,                 &      ! Number of subdomains
       nmate,                 &      ! # of materials
       nfiel,                 &      ! Number of fields
       kfl_field(7,mfiel),    &      ! Fields dimensions
       mcodb,                 &      ! Max # boundary codes
       ncodb                         ! # boundary codes
  
#ifdef NDIMEPAR
  ! much more comfortable this way; you can have a forder unix2d with -DTWODIM in the config.in
#ifdef TWODIM
  integer(ip), parameter   :: ndime = 2  ! # of space dimensions
#else
  integer(ip), parameter   :: ndime = 3  ! # of space dimensions
#endif
#else
  integer(ip)              :: ndime      ! # of space dimensions
#endif
  
  !------------------------------------------------------------------------
  !
  ! Parallelization
  !
  !------------------------------------------------------------------------
  
  integer(ip)              :: &
       npoi1,                 &      ! Internal nodes end node
       npoi2,                 &      ! Own-boundary start node
       npoi3,                 &      ! Own-boundary end node
       nedg1,                 &      ! Internal nodes end edge
       nedg2,                 &      ! Own-boundary start edge
       nedg3,                 &      ! Own-boundary end edge
       nelem_2,               &      ! nelem + fringe elements
       nboun_2,               &      ! nboun + fringe boundaries
       npoin_2,               &      ! npoin + fringe nodes
       npoin_own,             &      ! Own nodes
       npoin_halo,            &      ! Number of nodes up to halo nodes
       npoin_origi,           &      ! Number total of points (boundary not replicated)
       npoin_total,           &      ! Number total of points (boundary replicated)
       nelem_total,           &      ! Number total of elements
       nboun_total                   ! Number total of boundaries

end module def_kintyp_dims
  
