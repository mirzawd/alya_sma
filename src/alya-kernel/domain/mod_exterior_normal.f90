!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    extnor.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Compute exterior normals
!> @details Compute exterior normals to the computaitonal domain
!>          end identify boundary nodes. Internal boundary nodes,
!>          prescribed in the boundary connectivity field are also
!>          taken into account.
!>
!> @}
!-----------------------------------------------------------------------

module mod_exterior_normal

  ! Input
  use def_kintyp,            only : ip,rp,lg
  use def_elmtyp,            only : BAR3D,SHELL
  use def_elmtyp,            only : ELEXT
  use def_elmtyp,            only : BOEXT,BOINT
  use def_master,            only : INOTEMPTY
  use def_master,            only : INOTMASTER
  use def_master,            only : INOTSLAVE
  use def_master,            only : ISLAVE,intost
  use def_master,            only : NPOIN_TYPE
  use def_master,            only : lninv_loc,zeror
  use def_domain,            only : ndime,nboun
  use def_domain,            only : npoin,ltypb,npoin_own
  use def_domain,            only : nelem,nboun_2
  use def_domain,            only : vmass,memor_dom
  use def_domain,            only : nnode,lnodb
  use def_domain,            only : lelch,ngaus
  use def_domain,            only : ltype,elmar
  use def_domain,            only : lnods,coord
  use def_domain,            only : mnode,mnodb
  use def_domain,            only : lelbo,lnnob
  use def_domain,            only : ndimb,lnnod
  use def_kermod,            only : ndivi
  use mod_domain,            only : domain_memory_allocate
  use mod_domain,            only : domain_memory_deallocate
  use mod_maths,             only : maths_local_orthonormal_basis
  use def_coupli,            only : mcoup
  use mod_communications,    only : PAR_SUM
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_GHOST_BOUNDARY_EXCHANGE
  ! Output
  use def_domain,            only : nbopo,lpoty
  use def_domain,            only : lboch,exnor
  use def_domain,            only : meshe
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_messages,          only : messages_live
  use mod_optional_argument, only : optional_argument
  use mod_strings,           only : integer_to_string
  use mod_std
  use mod_bouder
  
  implicit none
  private

  public :: extnor
  public :: exterior_normal_destructor

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-01
  !> @brief   Destructor
  !> @details Exterior normal destructor
  !> 
  !-----------------------------------------------------------------------
  
  subroutine exterior_normal_destructor()

    call domain_memory_deallocate('EXNOR')
    
  end subroutine exterior_normal_destructor
  
  !-----------------------------------------------------------------------
  !
  !> @date    17/10/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Compute exterior normal
  !> @details Compute exterior normal EXNOR and determine
  !>          boundary nodes in LPOTY.
  !>
  !>          \verbatim
  !>
  !>          ITASK = 1 ... EXNOR and LPOTY
  !>                = 2 ... EXNOR
  !>
  !>          MESHE structure is also updated.
  !>
  !>          Let u = 1 on node, 0 elsewehere
  !>          Let We be the patch of node i
  !>          Let Se be the boundary of the patch We
  !>          Let n_x^i be the normal in x for node i
  !>
  !>          Interior node    Boundary node
  !>          0----0----0      0----1----0
  !>          |    |    |      |    |    |
  !>          0----1----0      0----0----0
  !>          |    |    |      |    |    |
  !>          0----0----0      0----0----0
  !>
  !>          n_x^i = Sum_e Sum_g dNi/dx |J|w_g = int_We du/dx = int_Se u.nx
  !>          So n_x^i  = 0 for interior node
  !>                   /= 0 for boundary nodes
  !>
  !>          Si definimos el valor nodal j de la componente i de la normal como:
  !>
  !>          ni(j) = \int_V dNj/dxi dv = \int_S N(j) ni ds
  !>
  !>          Entonces si calculamos mass = \sum_j ui(j).ni(j), tenemos
  !>
  !>          mass = \sum_j ui(j) . \sum_j ui(j).ni(j) = \int_S \sum_j Ui(j)N(j).ni ds
  !>                                                = \int_S u.n ds
  !>          Por lo tanto: \int_S u.n ds = 0 <=> mass = 0
  !>
  !>          OUTPUT
  !>             NBOPO ...................... Number of points located on the boundary
  !>             EXNOR(NDIME,NDIME,NBOPO) ... Local basis at boundary nodes
  !>             LPOTY(NPOIN) ............... defined by LPOTY(IPOIN):
  !>                                          = 0      if IPOIN is interior
  !>                                          = IBOPO  # of boundary point
  !>
  !>          Procedure:
  !>          1. Look for boundary nodes using EXWOR computed in entire domain:
  !>             - EXWOR == 0 ... Interior node
  !>             - EXWOR >  0 ... Boundary node
  !>          2. Loop over boundaries to detect possible boundary nodes (internal)
  !>             identified in first step as internal. For these nodes, we compute
  !>             ni(j) = \int_S N(j) ni ds
  !>          3. Still EXWOR can be null even for these nodes!
  !>          \endverbatim
  !
  !-----------------------------------------------------------------------

  subroutine extnor(itask)

    integer(ip), intent(in) :: itask
    integer(ip)             :: pblty,ipoin,iboun,inodb
    integer(ip)             :: ierro1,ierro2
    integer(ip)             :: kfl_internal
    integer(ip)             :: kfl_trailing_edge
    integer(ip)             :: dummi(4)
    real(rp)                :: direction(3)

    integer(ip), pointer    :: boundary_mask(:)
    real(rp),    pointer    :: exaux(:,:)
    real(rp),    pointer    :: exaux_int(:,:)

    call messages_live('COMPUTE EXTERIOR NORMALS')

    ierro1            = 0
    ierro2            = 0
    dummi             = 0
    kfl_internal      = 0
    kfl_trailing_edge = 0
    nullify(exaux)
    nullify(exaux_int)
    nullify(boundary_mask)

    if( INOTEMPTY ) then
       !
       ! Memory allocation and initializations
       !
       if( itask == 1 ) then
          call memory_alloca(memor_dom,'LPOTY','extnor',lpoty,npoin)
       end if
       call memory_alloca(memor_dom,'EXAUX','extnor',exaux,ndime,npoin)

       !-------------------------------------------------------------------
       !
       ! Compute element integrals that contribute to exaux
       !
       !-------------------------------------------------------------------

       call extnor_element_operations(exaux)

       !-------------------------------------------------------------------
       !
       ! NBOPO: Counting of the boundary points
       !
       !-------------------------------------------------------------------

       if( itask == 1 ) then

          nbopo = 0

          do ipoin = 1,npoin
             if( extnor_boundary_node(ndime,exaux(:,ipoin),vmass(ipoin),ipoin)) then
                nbopo        = nbopo + 1
                lpoty(ipoin) = 1
             end if
          end do

          !-------------------------------------------------------------------
          !
          ! Look for additional boundary nodes not detected by exaux.
          ! They can be given explicitly by the user in LNODB.
          ! The nodes of hole boundaries are taken off from this list.
          ! This is useful for example for Nastin which wants to automatically
          ! prescribe the pressure when LPOTY(IPOIN) /= 0.
          !
          ! This can be erroneous when boundaries have been explicitly declared
          ! on periodic or internal faces. Then Alya may not find local basis
          ! in setext.
          !
          !
          ! Mark these boundaries as internal boundaries
          !
          !-------------------------------------------------------------------

          do iboun = 1,nboun_2
             pblty = ltypb(iboun)
             if( lboch(iboun) /= BOEXT ) then
                do inodb = 1,nnode(pblty)
                   ipoin = lnodb(inodb,iboun)
                   if( ipoin <= npoin ) then
                      if( lpoty(ipoin) <= 0 ) then
                         if( lpoty(ipoin) == 0 ) then
                            nbopo        =  nbopo + 1
                            lpoty(ipoin) = -1
                         end if
                         lboch(iboun) = BOINT
                         kfl_internal = 1
                      end if
                   end if
                end do
             end if
          end do

       end if

    end if

    !-------------------------------------------------------------------
    !
    ! Compute normal for internal boundaries
    !
    !-------------------------------------------------------------------

    call PAR_MAX(kfl_internal)

    if( kfl_internal > 0 .and. INOTMASTER ) then

       call PAR_GHOST_BOUNDARY_EXCHANGE(lboch,'SUBSTITUTE','IN MY CODE')
       call memory_alloca(memor_dom,'EXAUX_INT',    'extnor',exaux_int,ndime,npoin)
       call memory_alloca(memor_dom,'BOUNDARY_MASK','extnor',boundary_mask,nboun)
       !
       ! Mark internal boundaries
       !
       do iboun = 1,nboun
          if( lboch(iboun) == BOINT ) then
             boundary_mask(iboun) = 1
          end if
       end do
       !
       ! Compute EXNOR on nodes belonging to internal boundaries
       !
       call extnor_boundary_operations(exaux_int,boundary_mask) !,'DO NOT IMPOSE PERIODICITY')
       !
       ! Subsititute exterior normal on nodes not already marked as boundary nodes
       !
       do ipoin = 1,npoin
          if( lpoty(ipoin) == -1 ) then
             exaux(1:ndime,ipoin) = exaux_int(1:ndime,ipoin)
          end if
       end do

       call memory_deallo(memor_dom,'EXAUX_INT',    'extnor',exaux_int)
       call memory_deallo(memor_dom,'BOUNDARY_MASK','extnor',boundary_mask)

    end if

    !-------------------------------------------------------------------
    !
    ! Now we should be only left with this kind of case:
    !
    ! o-----------o-----------o
    ! |           |           |
    ! |           |     i     |
    ! |           |___________o
    ! o-----------*___________o
    ! |           |     j     |
    ! |           |           |
    ! |           |           |
    ! o-----------o-----------o
    !
    ! In this case, the normal at node * would be only computed using
    ! i or j. Thus we should select only one side! How can we do this?
    !
    !-------------------------------------------------------------------

    if( INOTMASTER ) then
       do ipoin = 1,npoin
          if( lpoty(ipoin) == -1 ) then
             if( extnor_interior_node(ndime,exaux(:,ipoin),vmass(ipoin)) ) then
                lpoty(ipoin)      = -2
                kfl_trailing_edge =  1
             end if
          end if
       end do
    end if

    call PAR_MAX(kfl_trailing_edge)

    if( kfl_trailing_edge > 0 .and. INOTMASTER ) then

       call memory_alloca(memor_dom,'EXAUX_INT',    'extnor',exaux_int,ndime,npoin)
       call memory_alloca(memor_dom,'BOUNDARY_MASK','extnor',boundary_mask,nboun)
       !
       ! Mark trailing edge boundaries
       !
       do iboun = 1,nboun
          if( lboch(iboun) /= BOEXT ) then
             inodb_loop: do inodb = 1,lnnob(iboun)
                ipoin = lnodb(inodb,iboun)
                if( lpoty(ipoin) == -2 ) then
                   boundary_mask(iboun) = 1
                   exit inodb_loop
                end if
             end do inodb_loop
          end if
       end do
       !
       ! Compute EXNOR on nodes belonging to internal boundaries
       !
       direction = 0.72_rp
       call extnor_boundary_operations(exaux_int,boundary_mask,DIRECTION=direction)
       !
       ! Subsititute exterior normal on nodes not already marked as boundary nodes
       !
       do ipoin = 1,npoin
          if( lpoty(ipoin) == -2 ) then
             exaux(1:ndime,ipoin) = exaux_int(1:ndime,ipoin)
          end if
       end do

       call memory_deallo(memor_dom,'EXAUX_INT',    'extnor',exaux_int)
       call memory_deallo(memor_dom,'BOUNDARY_MASK','extnor',boundary_mask)

    end if

    !-------------------------------------------------------------------
    !
    ! Count internal and trailing edge nodes
    !
    !-------------------------------------------------------------------

    if( INOTEMPTY ) then
       dummi(1) = count(lpoty(1:npoin_own)==-1,KIND=ip)
       dummi(2) = count(lpoty(1:npoin_own)==-2,KIND=ip)
    end if

    !-------------------------------------------------------------------
    !
    ! Compute the vectors EXNOR (local basis) and LPOTY
    !
    !-------------------------------------------------------------------

    if( nbopo > 0 .and. INOTMASTER ) then
       if( itask == 1 ) call domain_memory_allocate('EXNOR')
       call extnor_normalize_exterior_normal(itask,exaux,ierro1,ierro2)
    end if
    if( INOTMASTER ) call memory_deallo(memor_dom,'EXAUX','extnor',exaux)
    !-------------------------------------------------------------------
    !
    ! Check errors
    !
    !-------------------------------------------------------------------

    dummi(3) = ierro1
    dummi(4) = ierro2
    call PAR_SUM(4_ip,dummi)
    kfl_internal      = dummi(1)
    kfl_trailing_edge = dummi(2)
    ierro1            = dummi(3)
    ierro2            = dummi(4)

    if( ierro1            > 0 ) call messages_live('SOME NORMALS ARE NOT DEFINED: THEY MAY BE ON A COUPLING INTERFACE','WARNING')
    if( ierro2            > 0 ) call messages_live('SOME LOCAL BASIS ARE NOT DEFINED','WARNING')
    if( kfl_internal      > 0 ) call messages_live('EXTERIOR NORMAL: NODES ON INTERNAL BOUNDARIES HAVE BEEN FOUND= '//trim(intost(kfl_internal)),'WARNING')
    if( kfl_trailing_edge > 0 ) call messages_live('EXTERIOR NORMAL: TRAILING EDGE NODES HAVE BEEN FOUND= '//trim(intost(kfl_trailing_edge)),'WARNING')

    !-------------------------------------------------------------------
    !
    ! Fill in mesh type
    !
    !-------------------------------------------------------------------

    meshe(ndivi) % exnor => exnor
    meshe(ndivi) % nbopo =  nbopo
    meshe(ndivi) % lpoty => lpoty


  end subroutine extnor

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/11/2017
  !> @brief   Compute exterior normals
  !> @details Computes exterior normals EXAUX(NDIME,NPOIN)
  !>
  !-----------------------------------------------------------------------

  subroutine extnor_element_operations(exaux)

    real(rp),    intent(inout), pointer           :: exaux(:,:)        !< Exterior normal
    integer(ip)                                   :: ielem,igaus,pelty
    integer(ip)                                   :: inode,ifoun,ipoin
    integer(ip)                                   :: knode,idime,pnode
    real(rp)                                      :: gpcar(ndime,mnode)
    real(rp)                                      :: xjaci(ndime,ndime)
    real(rp)                                      :: xjacm(ndime,ndime)
    real(rp)                                      :: elcod(ndime,mnode)
    real(rp)                                      :: elaux(ndime,mnode)
    real(rp)                                      :: gpvol,gpdet

    if( INOTEMPTY ) then

       !-------------------------------------------------------------------
       !
       ! Compute element integrals that contribute to exaux
       !
       !-------------------------------------------------------------------
       !
       !$OMP  PARALLEL DO SCHEDULE (STATIC)                           &
       !$OMP  DEFAULT (NONE)                                          &
       !$OMP  PRIVATE ( ielem, pelty, pnode, inode, ipoin, igaus,     &
       !$OMP            gpcar, gpdet, xjacm, xjaci, gpvol, elcod,     &
       !$OMP            elaux, ifoun, knode )                         &
       !$OMP   SHARED ( ltype, nnode, lnods, coord, nelem, ngaus,     &
       !$OMP            lelch, elmar,                                 &
#ifndef NDIMEPAR
       !$OMP            ndime,                                        &
#endif
       !$OMP            exaux )
       !
       do ielem = 1,nelem
          pelty = ltype(ielem)

          ifoun = 1

          if( pelty > 0 .and. ifoun == 1 ) then
             pnode = nnode(pelty)
             if( lelch(ielem) == ELEXT ) then
                knode = 1
             else
                knode = pnode
             end if
             elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
             elaux = 0.0_rp
             if( pelty == BAR3D ) then
                elaux(:,1) = elcod(:,1)-elcod(:,2)
                elaux(:,2) = elcod(:,2)-elcod(:,1)
             else
                do igaus = 1,ngaus(pelty)
                   call jacobi(&
                        ndime,pnode,elcod,elmar(pelty) % deriv(1,1,igaus),&
                        xjacm,xjaci,gpcar,gpdet)
                   gpvol = elmar(pelty) % weigp(igaus) * gpdet
                   elaux(1:ndime,1:knode) = elaux(1:ndime,1:knode) + gpcar(1:ndime,1:knode) * gpvol
                end do
             end if
             do inode = 1,knode
                ipoin = lnods(inode,ielem)
                do idime = 1,ndime
                   !$OMP ATOMIC
                   exaux(idime,ipoin) = exaux(idime,ipoin) + elaux(idime,inode)
                end do
             end do
          end if

       end do

       !-------------------------------------------------------------------
       !
       ! Modify EXAUX due to periodicity and Parall service
       !
       !-------------------------------------------------------------------

       call rhsmod(ndime,exaux)

    end if

  end subroutine extnor_element_operations

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/11/2017
  !> @brief   Compute exterior normals
  !> @details Computes exterior normals EXAUX(NDIME,NPOIN)
  !>
  !-----------------------------------------------------------------------

 subroutine extnor_boundary_operations(exaux,boundary_mask,what,DIRECTION)

    real(rp),    intent(inout), pointer           :: exaux(:,:)         !< Exterior normal
    integer(ip), intent(in),    pointer, optional :: boundary_mask(:)   !< Boundary mask
    character(*),intent(in),             optional :: what
    real(rp),    intent(in),             optional :: DIRECTION(3)
    integer(ip)                                   :: igaub,pnodb,pblty
    integer(ip)                                   :: ielem,pnode,ipoin
    integer(ip)                                   :: inodb,iboun
    logical(lg)                                   :: compute_boundary
    logical(lg)                                   :: impose_periodicity
    real(rp)                                      :: gbsur,eucta,dotpro
    real(rp)                                      :: bocod(ndime,mnodb)
    real(rp)                                      :: elcod(ndime,mnode)
    real(rp)                                      :: baloc(ndime,ndime)
    real(rp)                                      :: boaux(ndime,mnodb)

    if( INOTMASTER ) then

       do iboun = 1,nboun

          compute_boundary = .true.
          if( present(boundary_mask) ) then
             if( boundary_mask(iboun) == 0 ) compute_boundary = .false.
          end if

          if( compute_boundary ) then

             pnodb = lnnob(iboun)
             pblty = abs(ltypb(iboun))
             ielem = lelbo(iboun)
             pnode = lnnod(ielem)
             boaux = 0.0_rp

             bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun))
             elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))

             gauss_points: do igaub = 1,ngaus(pblty)
                call bouder(&
                     pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),&
                     bocod,baloc,eucta)
                gbsur = elmar(pblty) % weigp(igaub)*eucta
                call chenor(pnode,baloc,bocod,elcod)                      ! Check normal
                if( present(DIRECTION) ) then
                   dotpro = dot_product(direction(1:ndime),baloc(1:ndime,ndime))
                   if( dotpro < 0.0_rp ) baloc = 0.0_rp
                end if
                do inodb = 1,pnodb
                   boaux(1:ndime,inodb) = boaux(1:ndime,inodb) &
                        + elmar(pblty) % shape(inodb,igaub) * gbsur * baloc(1:ndime,ndime)
                end do

             end do gauss_points

             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                exaux(1:ndime,ipoin) = exaux(1:ndime,ipoin) + boaux(1:ndime,inodb)
             end do

          end if

       end do

       !-------------------------------------------------------------------
       !
       ! Modify EXAUX due to periodicity and Parall service
       !
       !-------------------------------------------------------------------
       !
       ! Option 'DO NOT IMPOSE PERIODICITY' is usefull not to have zero normal
       ! on periodic nodes
       !
       impose_periodicity = .true.
       if( present(what) ) then
          if( trim(what) == 'DO NOT IMPOSE PERIODICITY') then
             impose_periodicity = .false.
          end if
       end if

       if( impose_periodicity ) then
          call rhsmod(ndime,exaux)
       else if( INOTEMPTY ) then
          call pararr('SLX',NPOIN_TYPE,npoin*ndime,exaux)
       end if

    end if

  end subroutine extnor_boundary_operations

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/11/2017
  !> @brief   Normalize exterior normal and define boundary nodes
  !> @details Normalize exterior normal and define boundary nodes:
  !>          If ITASK = 1, compute both EXNOR and LPOTY.
  !>          If ITASK = 2, then update only EXNOR. This is called when
  !>          the mesh has changed for example due to a mesh movement
  !>          by Alefor.
  !>          \verbatim
  !>          EXNOR(NDIME,NDIME,NBOPO) ... Local basis at boundary nodes
  !>          LPOTY(NPOIN) ............... defined by LPOTY(IPOIN)
  !>                                       = 0     If IPOIN is interior
  !>                                       = IBOPO Number of boundary point
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine extnor_normalize_exterior_normal(itask,exwor,ierro1,ierro2)

    !   To debug only ->
    use mod_output,         only : output_mesh_gid_format
    use def_kermod,         only : ndivi
    use def_master,         only : lun_outpu
    ! <- o debug only

    integer(ip), intent(in)            :: itask       !< Task
    real(rp),    intent(in),  pointer  :: exwor(:,:)  !< Exterior normal on all nodes
    integer(ip), intent(out), optional :: ierro1      !< Error message for normals
    integer(ip), intent(out), optional :: ierro2      !< Error message tangent basis
    integer(ip)                        :: ipoin,ibopo
    integer(ip)                        :: jerro,ii
    real(rp)                           :: xnorm
    logical(lg)                        :: debug_mode
    logical(lg),              pointer  :: mark_npoin(:)

    debug_mode = .false.

    if( present(ierro1) ) ierro1 = 0
    if( present(ierro2) ) ierro2 = 0

    select case ( itask )

    case ( 1_ip )

       !-------------------------------------------------------------------
       !
       ! EXNOR and LPOTY
       !
       !-------------------------------------------------------------------

       exnor = 0.0_rp
       ibopo = 0
       do ipoin = 1,npoin
          if( abs(lpoty(ipoin)) >= 1 ) then
             xnorm = sqrt(dot_product(exwor(1:ndime,ipoin),exwor(1:ndime,ipoin)))
             if( xnorm > 0.0_rp ) then
             !if( extnor_boundary_node(ndime,exwor(:,ipoin),vmass(ipoin),ipoin) ) then
                xnorm                  = sqrt(dot_product(exwor(1:ndime,ipoin),exwor(1:ndime,ipoin)))
                ibopo                  = ibopo+1
                lpoty(ipoin)           = ibopo
                exnor(1:ndime,1,ibopo) = exwor(1:ndime,ipoin)/xnorm ! Normalize normal
             else
                lpoty(ipoin)           = -1
             end if
          end if
       end do
       !
       ! DEBUG_MODE: output mesh patch around boundary nodes without a normal
       !
       if( ibopo /= nbopo .and. debug_mode ) then
          nullify(mark_npoin)
          call memory_alloca(memor_dom,'MARK_NPOIN','extnor',mark_npoin,npoin)
          do ipoin = 1,npoin
             if( lpoty(ipoin) == -1 ) mark_npoin(ipoin) = .true.
          end do
          call output_mesh_gid_format(meshe(ndivi),'EXNOR_PROBLEM',lun_outpu,mark_npoin,RESULT_INT=lninv_loc)
          call memory_deallo(memor_dom,'MARK_NPOIN','extnor',mark_npoin)
       end if
       !
       ! Look for additional boundary nodes not detected by exwor and
       ! explicitely declared by lnodb: These are for example wet nodes
       ! of coupling. But normally, the exterior normal is not required
       ! anywhere for this nodes, so we put any value... just in case.
       !
       if( ibopo /= nbopo ) then
          do ipoin = 1,npoin
             if( lpoty(ipoin) <= -1 ) then
                ibopo        = ibopo + 1
                lpoty(ipoin) = ibopo
                if( present(ierro1) ) ierro1 = ierro1 + 1
                exnor(:,1,ibopo) = 0.0_rp
                exnor(1,1,ibopo) = 1.0_rp
                if( mcoup == 0 .or. debug_mode ) then
                   ii = lninv_loc(ipoin)
                   call messages_live('SETEXT: NORMAL COULD NOT BE COMPUTED FOR NODE '//integer_to_string(ii)//'; ARBITRARILY SET TO 1,0,0 ','WARNING')
                end if
             end if
          end do
       end if

    case ( 2_ip )

       !-------------------------------------------------------------------
       !
       ! EXNOR only
       !
       !-------------------------------------------------------------------

       exnor = 0.0_rp
       do ipoin = 1,npoin
          if( lpoty(ipoin) /= 0 ) then
             ibopo = lpoty(ipoin)
             xnorm = dot_product(exwor(1:ndime,ipoin),exwor(1:ndime,ipoin))
             if( xnorm > 0.0_rp ) then
                exnor(1:ndime,1,ibopo) = exwor(1:ndime,ipoin)/sqrt(xnorm) ! Normalize normal
             else
                if( present(ierro1) ) ierro1 = ierro1 + 1
                exnor(:,1,ibopo) = 0.0_rp
                exnor(1,1,ibopo) = 1.0_rp
             end if
          end if
       end do

    end select

    !--------------------------------------------------------------------
    !
    ! Construct the tangent vectors to define the local basis
    !
    !--------------------------------------------------------------------

    do ipoin = 1,npoin
       ibopo = lpoty(ipoin)
       if( ibopo > 0 ) then
          jerro = 0
          call maths_local_orthonormal_basis(ndime,exnor(:,:,ibopo),jerro)
          if( present(ierro2) ) ierro2 = ierro2 + jerro
          if( jerro /= 0 ) then
             write(*,*) 'SETEXT: TANGENT COULD NOT BE COMPUTED FOR NODE ',lninv_loc(ipoin)
          end if
       end if
    end do

  end subroutine extnor_normalize_exterior_normal


  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/11/2017
  !> @brief   Check if a node is interior
  !> @details Check if a node is interior
  !>
  !-----------------------------------------------------------------------

  function extnor_interior_node(kdime,exaux,vmass)

    integer(ip), intent(in) :: kdime
    real(rp),    intent(in) :: exaux(kdime)
    real(rp),    intent(in) :: vmass
    logical(lg)             :: extnor_interior_node

    extnor_interior_node = .not. extnor_boundary_node(kdime,exaux,vmass)

  end function extnor_interior_node

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/11/2017
  !> @brief   Check if a node is on boundary
  !> @details Check if a node is a boundary or interior node
  !>
  !-----------------------------------------------------------------------

  function extnor_boundary_node(kdime,exaux,vmass,IPOIN,TOLER)

    integer(ip), intent(in)           :: kdime                !< Dimension
    real(rp),    intent(in)           :: exaux(kdime)         !< Non-normalized normal
    real(rp),    intent(in)           :: vmass                !< Mass
    integer(ip), intent(in), optional :: IPOIN                !< Node, for debug purpose
    real(rp),    intent(in), optional :: TOLER                !< Tolerance
    logical(lg)                       :: extnor_boundary_node
    real(rp)                          :: exwn2,my_toler

    my_toler = optional_argument(1.0e-6_rp,TOLER)
    extnor_boundary_node = .false.

    select case ( kdime )

    case ( 3_ip )

       exwn2 = sqrt(dot_product(exaux(1:3),exaux(1:3)))
       if( exwn2 > my_toler*vmass**(1.0_rp/3.0_rp) ) extnor_boundary_node = .true.

    case ( 2_ip )

       exwn2 = sqrt(dot_product(exaux(1:2),exaux(1:2)))
       if( exwn2 > my_toler*sqrt(vmass) ) extnor_boundary_node = .true.

    case ( 1_ip )

       exwn2 = exaux(1) * exaux(1)
       if( exwn2 > my_toler ) extnor_boundary_node = .true.

    case default

       call runend('EXTNOR_BOUNDARY_NODE: WRONG DIMENSION')

    end select

  end function extnor_boundary_node

end module mod_exterior_normal
!> @}
