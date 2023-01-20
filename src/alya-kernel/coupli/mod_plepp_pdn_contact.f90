!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    mod_plepp_pdn_contact.f90
!> @author  gguillam
!> @date    2022-05-25
!> @brief   PLE++ toolbox for PDN contact algorithm
!> @details 
!>          Used by Solidz module
!>          This is a clean version from M.Zavala modules particularly
!>          focused on the use of the PDN contact algorithm using PLE++
!>
!>          Alya use the following public variables/subroutines:
!>
!>          PLEPP_CPLNG
!>          CNT_SENDRECV
!>          CNT_CPLNG                     /src/alya-kernel/domain/domain.f90
!>                                        /src/alya-core/mod_coupling_driver.f90
!>          plepp_pdn_init                /src/alya-kernel/parall/par_code_split_universe.f90
!>          plepp_pdn_driver_init         /src/alya-kernel/domain/domain.f90
!>          plepp_pdn_driver_memall       /src/alya-kernel/domain/domain.f90
!>          plepp_pdn_driver_sendrecv     /src/alya-core/mod_coupling_driver.f90
!>          plepp_pdn_driver_exchange02   /src/modules/solidz/mod_sld_pdn_contact.f90       
!>          plepp_pdn_dynamic_reduce_sum  /src/modules/solidz/mod_sld_pdn_contact.f90
!>
!> @{
!-----------------------------------------------------------------------

module mod_plepp_pdn_contact

#if defined COMMDOM && COMMDOM == 2

  use def_kintyp, only : ip,rp,lg
  use def_master, only : INOTMASTER,IMASTER,ISEQUEN, islave, inotslave, iparall
  use def_master, only : routp, dtime, kfl_paral
  use def_master, only : ittim, cutim, mitim
  use def_kermod, only : kfl_conta
  use def_domain, only : mnode, nelem, ndime, npoin, nboun
  use def_domain, only : kfl_codno
  use def_domain, only : ltype, nnode, ngaus, lnods, coord, lpoty, lelbo 
  use mod_parall, only : PAR_COMM_CURRENT,PAR_UNIVERSE_SIZE
  use mod_parall, only : PAR_COMM_WORLD,PAR_COMM_UNIVERSE,PAR_COMM_MY_CODE

  implicit none

  logical(lg)              :: initiated = .false.
  integer(ip), parameter   :: n_times = 12_ip
  logical(lg)              :: CNT_SENDRECV(n_times) = .false.
  logical(lg)              :: resid = .false.
  
  type COMMDOM_COUPLING
     integer(ip)           :: code_i        = -1_ip !< CODE.   coupling%code_source
     integer(ip)           :: code_j        = -1_ip
     integer(ip)           :: module_i      = -1_ip !< MODULE. coupling%module_sourc
     integer(ip)           :: module_j      = -1_ip
     integer(ip)           :: what_i        = -1_ip
     integer(ip)           :: what_j        = -1_ip
     logical(lg)           :: current_what(5_ip) = .false.
     integer(ip)           :: current_code  = -1_ip
     integer(ip)           :: n_dof         = -1_ip
     integer(ip)           :: n_pts         = -1_ip
     real(rp),   pointer   :: var_ji(:,:)
     real(rp),   pointer   :: var_ij(:,:)
  end type COMMDOM_COUPLING

  type(COMMDOM_COUPLING) CNT_CPLNG
  
  type COMMDOM_PLEPP_COUPLING
     integer(ip)           :: local_comm           = -1_ip
     integer(ip)           :: commij               = -1_ip
     integer(ip)           :: commji               = -1_ip
     integer(ip)           :: n_ij                 = -1_ip
     integer(ip)           :: n_ji                 = -1_ip
     integer(ip)           :: stride               = -1_ip
     real(rp)              :: tol                  =  0.0_rp
     real(rp),    pointer  :: var_ij(:,:)          => null()
     real(rp),    pointer  :: var_ji(:,:)          => null()
     real(rp),    pointer  :: tetra_coords_j(:,:)  => null()
     real(rp),    pointer  :: dist_coords_j(:)     => null()
     integer(ip), pointer  :: dist_locations_i(:)  => null()
     integer(ip), pointer  :: interior_list_j(:)   => null()
     real(rp),    pointer  :: dist_props_j(:)      => null()
     character(64)         :: namei       = ''
     character(64)         :: namej       = ''
     character(64)         :: app_type    = ''
     character(64)         :: app_name    = ''
     character(64)         :: module_name = ''
  end type COMMDOM_PLEPP_COUPLING

  type(COMMDOM_PLEPP_COUPLING) PLEPP_CPLNG

  type COMMDOM_VERTEX_PROPERTIES
     integer(ip)           :: ndime_props
     integer(ip)           :: npoin    
     integer(ip)           :: nelem 
     integer(ip), pointer  :: ltype(:)   => null()
     integer(ip), pointer  :: lnods(:,:) => null()
     integer(ip), pointer  :: lpoin(:)   => null()
     real(rp),    pointer  :: coord(:,:) => null()
     real(rp),    pointer  :: displ(:,:) => null() 
  end type COMMDOM_VERTEX_PROPERTIES

  type(COMMDOM_VERTEX_PROPERTIES) CPLNG_PROPS

  public :: PLEPP_CPLNG
  public :: CNT_SENDRECV
  public :: CNT_CPLNG
  public :: plepp_pdn_init                
  public :: plepp_pdn_driver_init         
  public :: plepp_pdn_driver_memall       
  public :: plepp_pdn_driver_sendrecv     
  public :: plepp_pdn_driver_exchange02        
  public :: plepp_pdn_dynamic_reduce_sum  
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_init()
    use def_master, only : intost
    implicit none

    integer(ip)    :: world_comm
    character(64)  :: token
    integer(ip)    :: oki, okj
    integer(4)     :: iarg

    PLEPP_CPLNG%module_name =  ''    
    PLEPP_CPLNG%app_type = "SYRTHES 4" !< coupling with saturne

    PLEPP_CPLNG%namei = "FLUID"
    PLEPP_CPLNG%namej = "SOLID"

    PLEPP_CPLNG%tol   =  1.0e-3_rp

    call commdom_create()

    PLEPP_CPLNG%app_name = ""
    do iarg = 1, command_argument_count() !iargc()
       token = ""
       call get_command_argument(iarg, token)
       call commdom_set_argvs(trim(token), len_trim(token))
       if( token == '--name' ) then
          call get_command_argument(iarg+1, token)
          PLEPP_CPLNG%app_name = trim(token)
       end if
    enddo

    !For some reason, this parser from plepp does not work, commdom_get_argvs does not return the string
    !token = "--name"
    !call commdom_analyse_argvs(trim(token), len_trim(token))

    !PLEPP_CPLNG%app_name = ""
    !call commdom_get_argvs(PLEPP_CPLNG%app_name)
    call commdom_set_names(trim(PLEPP_CPLNG%app_type), len_trim(PLEPP_CPLNG%app_type), &
         trim(PLEPP_CPLNG%app_name), len_trim(PLEPP_CPLNG%app_name))
    PLEPP_CPLNG%module_name = trim(PLEPP_CPLNG%app_name)
    
    !-----------------------------------------------------------------------||---!

    !-----------------------------------------------------------------------||---!
    
#if defined USEMPIF08
    world_comm = PAR_COMM_UNIVERSE % MPI_VAL
#else
    world_comm = PAR_COMM_UNIVERSE 
#endif
    
    call commdom_create_commij(world_comm, PLEPP_CPLNG%local_comm)

#if defined USEMPIF08
    PAR_COMM_WORLD   % MPI_VAL = PLEPP_CPLNG%local_comm
    PAR_COMM_CURRENT % MPI_VAL = PLEPP_CPLNG%local_comm
#else
    PAR_COMM_WORLD             = PLEPP_CPLNG%local_comm
    PAR_COMM_CURRENT           = PLEPP_CPLNG%local_comm
#endif
    
    !-----------------------------------------------------------------------||---!

    !-----------------------------------------------------------------------||---!
    
    PLEPP_CPLNG%commij = -1
    PLEPP_CPLNG%commji = -1

    oki = -1
    okj = -1

    call commdom_who_iam(trim(PLEPP_CPLNG%namei), len_trim(PLEPP_CPLNG%namei), oki)
    if(oki==1) then
       call commdom_who_areyou(trim(PLEPP_CPLNG%namej), len_trim(PLEPP_CPLNG%namej), okj)
       if(okj==1) then
          call commdom_get_commij(trim(PLEPP_CPLNG%namej), len_trim(PLEPP_CPLNG%namej), PLEPP_CPLNG%commij)
       endif
    endif

    if( oki<0 .and. okj<0 ) then
       oki = -1
       okj = -1

       call commdom_who_iam(trim(PLEPP_CPLNG%namej), len_trim(PLEPP_CPLNG%namej), okj)
       if(okj==1) then
          call commdom_who_areyou(trim(PLEPP_CPLNG%namei), len_trim(PLEPP_CPLNG%namei), oki)
          if(oki==1) then
             call commdom_get_commij(trim(PLEPP_CPLNG%namei), len_trim(PLEPP_CPLNG%namei), PLEPP_CPLNG%commji)
          endif
       endif
    endif

    if( oki<0 .and. okj<0 ) then
       print *,"plepp_pdn_init: PLEPP failed to set up communication, oki="//trim(intost(oki))//", okj="//trim(intost(okj))
       stop 1
    end if


    !-----------------------------------------------------------------------||---!

    if( PLEPP_CPLNG%commij /= -1 .or. PLEPP_CPLNG%commji /= -1 ) then
       PAR_UNIVERSE_SIZE = -1 !=> if( PAR_UNIVERSE_SIZE > 1 ) ...
    endif
    
    !-----------------------------------------------------------------------||---!

  end subroutine plepp_pdn_init

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_compare_dtinv(dt_inv)

    implicit none

    real(rp), intent(inout) :: dt_inv
    real(rp)                :: send, recv, dtmin
    integer(ip)             :: comm, n_send, n_recv
    !!character(16)           :: str(3)

    if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then

       if( .not.(PLEPP_CPLNG%commij == -1) ) then
          comm = PLEPP_CPLNG%commij
       endif
       if( .not.(PLEPP_CPLNG%commji == -1) ) then
          comm = PLEPP_CPLNG%commji
       endif

       send   =  dt_inv
       recv   = -1.0_rp
       n_send =  1
       n_recv =  1
       call commdom_sendrecv_real(send, n_send, recv, n_recv, PLEPP_CPLNG%local_comm, comm)
       call commdom_bcast_real(recv, n_recv, PLEPP_CPLNG%local_comm, comm)

       dtmin  = min(1.0_rp/send, 1.0_rp/recv)
       dt_inv = 1.0_rp/dtmin

       !if( IMASTER .or. ISEQUEN ) then
       !   write(str(1),'(e11.4)') 1.0_rp/send
       !   write(str(2),'(e11.4)') 1.0_rp/recv
       !   write(str(3),'(e11.4)') dtmin
       !   print *, ""
       !   print *, '[commdom_plepp_compare_dtinv] '// trim(PLEPP_CPLNG%module_name) //' '&
       !        //'['//trim(str(1))&
       !        //','//trim(str(2))&
       !        //','//trim(str(3))&
       !        //']'
       !endif

    endif

  end subroutine plepp_pdn_compare_dtinv

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_locator_send_nodal_var00(&
       elemts_i, coords_j, n_dist_j, tetracoords_j)

    implicit none

    integer(ip), intent(in)  :: n_dist_j
    integer(ip), intent(in)  :: elemts_i(n_dist_j)
    real(rp),    intent(in)  :: coords_j(n_dist_j*ndime) 
    real(rp),    intent(out) :: tetracoords_j(n_dist_j,mnode)

    integer(ip)              :: ii, ielem, pelty, pnode
    integer(ip)              :: vertices_i(mnode)
    real(rp)                 :: vol_coords_j(mnode)
    real(rp)                 :: point_j(ndime)

    tetracoords_j(:,:) = 0.0_rp
    
    if( INOTMASTER ) then

       elementary_loop: do ii = 1,n_dist_j
          
          ielem               = elemts_i(ii)
          pelty               = ltype(ielem)
          pnode               = nnode(pelty)
          vertices_i(1:mnode) = -1
          vertices_i(1:pnode) = lnods(1:pnode,ielem)
          point_j(1:ndime)    = coords_j(ndime*ii-ndime+1:ndime*ii) 
          call commdom_locator_interpolation(&
               pelty, PLEPP_CPLNG%tol, coord(1:ndime,1:npoin), &
               vertices_i(1:pnode), point_j(1:ndime), vol_coords_j(1:pnode) )
          tetracoords_j(ii,1:pnode) = vol_coords_j(1:pnode)
          
       end do elementary_loop

    end if

  end subroutine plepp_pdn_locator_send_nodal_var00

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_locator_send_nodal_var01(&
       npoin_loc,props_i,elemts_i,n_dist_j,tetracoords_j, props_j)

    implicit none

    integer(ip), intent(in)  :: npoin_loc
    real(rp),    intent(in)  :: props_i(npoin_loc)
    integer(ip), intent(in)  :: n_dist_j
    integer(ip), intent(in)  ::      elemts_i(n_dist_j)
    real(rp),    intent(in)  :: tetracoords_j(n_dist_j,mnode)
    real(rp),    intent(out) ::       props_j(n_dist_j)

    integer(ip)              :: ii, ielem, pelty, pnode
    integer(ip)              :: vertices_i(mnode)
    real(rp)                 :: vol_coords_j(mnode)
    real(rp)                 :: prop_j(mnode)

    props_j(:) = 0.0_rp
    
    if( INOTMASTER ) then
       
       elementary_loop: do ii = 1,n_dist_j
          
          ielem                 = elemts_i(ii)
          pelty                 = ltype(ielem)
          pnode                 = nnode(pelty)
          vertices_i(  1:pnode) = lnods(1:pnode,ielem)
          vol_coords_j(1:pnode) = tetracoords_j(ii,1:pnode)
          prop_j(1:pnode)       = props_i( vertices_i(1:pnode) )
          props_j(ii)           = dot_product( prop_j(1:pnode), vol_coords_j(1:pnode) )
          
       enddo elementary_loop

    endif

  end subroutine plepp_pdn_locator_send_nodal_var01

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_dynamic_set_mesh(stride)

    implicit none

    integer(ip), optional, intent(in) :: stride

    if( .not.initiated ) then

       call plepp_pdn_dynamic_create(                  PLEPP_CPLNG        )
       call plepp_pdn_dynamic_allocate(                PLEPP_CPLNG,stride )
       call plepp_pdn_dynamic_get_coupling_properties( PLEPP_CPLNG        )
       
       initiated = .true.
       
    end if

  end subroutine plepp_pdn_dynamic_set_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_dynamic_mesh_size()

    implicit none
    
    integer(ip)          :: npoin_loc,kpoin,ipoin,iboun,ielem,pnode,inode
    integer(ip), pointer :: lpoin_aux(:) => null()

    if( INOTMASTER ) then

       allocate(lpoin_aux(npoin))
       lpoin_aux = 0_ip
       npoin_loc = 0_ip
       do iboun = 1,nboun
          ielem = lelbo(iboun)
          pnode = nnode(ltype(ielem))
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             do kpoin = 1,npoin
                if( ipoin /= lpoin_aux(kpoin) .and. lpoin_aux(kpoin) == 0_ip ) then
                   lpoin_aux(kpoin) = ipoin
                   npoin_loc = npoin_loc + 1
                end if
             end do
          end do
       end do

       CPLNG_PROPS%npoin = npoin_loc
       CPLNG_PROPS%nelem = nboun

       if( .not.associated(CPLNG_PROPS%coord) ) allocate(CPLNG_PROPS%coord(ndime,npoin_loc))
       if( .not.associated(CPLNG_PROPS%ltype) ) allocate(CPLNG_PROPS%ltype(nboun))
       if( .not.associated(CPLNG_PROPS%lpoin) ) allocate(CPLNG_PROPS%lpoin(npoin_loc))
       if( .not.associated(CPLNG_PROPS%lnods) ) allocate(CPLNG_PROPS%lnods(mnode,nboun))

       CPLNG_PROPS%coord = 0.0_rp
       CPLNG_PROPS%ltype = 0_ip
       CPLNG_PROPS%lnods = 0_ip
       CPLNG_PROPS%lpoin = 0_ip

       kpoin = 0
       do ipoin = 1,npoin
          if( lpoin_aux(ipoin) > 0 ) then
             kpoin = kpoin + 1
             CPLNG_PROPS%lpoin(kpoin) = ipoin 
             CPLNG_PROPS%coord(1:ndime,kpoin) = coord(1:ndime,ipoin)
          end if
       end do

       deallocate(lpoin_aux)

       do iboun = 1,nboun
          ielem = lelbo(iboun)
          pnode = nnode(ltype(ielem))
          CPLNG_PROPS%ltype(iboun)         = ltype(ielem)
          CPLNG_PROPS%lnods(1:pnode,iboun) = lnods(1:pnode,ielem)
       end do

    end if
    
  end subroutine plepp_pdn_dynamic_mesh_size
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_dynamic_create( PLEPP )

    use def_master,         only : displ
    use mod_communications, only : PAR_MAX

    implicit none

    type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP

    real(rp),    pointer :: vertex_coords_i(:,:) => null()
    integer(ip), pointer ::    vertex_num_i(:,:) => null()
    integer(ip), pointer ::   vertex_type_i(:  ) => null()
    real(rp),    pointer :: vertex_coords_j(:,:) => null()
    real(rp),    pointer :: vertex_props_j(:,:)  => null()

    integer(ip)          :: n_vertices_i, n_elements_i, n_vertices_j
    integer(ip)          :: element_type
    !
    ! Initialize
    !
    n_vertices_i =  0_ip
    n_elements_i =  0_ip
    n_vertices_j =  0_ip
    element_type = -1_ip

    if( PLEPP%commij /= -1 .or. PLEPP%commji /= -1 ) then
       !
       ! Init mesh
       !
       if( INOTMASTER ) then

          CPLNG_PROPS%ndime_props = ndime

          n_vertices_i    =  npoin
          n_elements_i    =  nelem          
          vertex_coords_i => coord
          vertex_num_i    => lnods
          vertex_type_i   => ltype

          n_vertices_j    =  npoin
          vertex_coords_j => coord
          vertex_props_j  => displ(:,:,1)

       else

          allocate( vertex_coords_i(0,0) )
          allocate(    vertex_num_i(0,0) )
          allocate(   vertex_type_i(0  ) )
          allocate( vertex_coords_j(0,0) )
          allocate(  vertex_props_j(0,0) )

          vertex_coords_i = 0_rp
          vertex_num_i    = 0_ip
          vertex_type_i   = 0_ip
          vertex_coords_j = 0.0_rp
          vertex_props_j  = 0.0_rp

       endif
       !
       ! Element type
       !
       if( all(vertex_type_i==30) .or. all(vertex_type_i==10) .or. all(vertex_type_i==37) ) then
          if( all(vertex_type_i==30) .or. all(vertex_type_i==10)) element_type = 4 ! TETRAS=30 | TRIA=10
          if( all(vertex_type_i==37)                            ) element_type = 6 ! HEXAS=37
       else
          call runend('COMMDOM: Hybrid mesh not ready')
       end if
       call PAR_MAX(element_type,'IN MY CODE')
       !
       ! Creat locator
       !
       if( .not.(PLEPP%commij == -1) ) then
          call commdom_locator_create2(PLEPP%local_comm, PLEPP%commij, PLEPP%tol, element_type)
       endif
       if( .not.(PLEPP%commji == -1) ) then
          call commdom_locator_create2(PLEPP%local_comm, PLEPP%commji, PLEPP%tol, element_type)
       endif

       call commdom_locator_set_cs_mesh(&
            n_vertices_i, &
            n_elements_i, &
            vertex_coords_i, &
            vertex_num_i, &
            vertex_type_i, &
            n_vertices_j, &
            vertex_coords_j, &
            ndime, &
            vertex_props_j, &
            CPLNG_PROPS%ndime_props)

    end if
    
  end subroutine plepp_pdn_dynamic_create

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_dynamic_allocate( PLEPP, stride )

    implicit none

    type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
    integer(ip),        optional, intent(in)    :: stride
    integer(ip)                                 :: n_recv, n_send

    if( PLEPP%commij /= -1 .or. PLEPP%commji /= -1  ) then

       n_recv = 0
       n_send = 0
       call commdom_locator_get_n_dist_points( n_send )
       call commdom_locator_get_n_interior(    n_recv )
       if( .not.associated( PLEPP%tetra_coords_j )  ) allocate( PLEPP%tetra_coords_j(            n_send, mnode) )
       if( .not.associated( PLEPP%dist_locations_i) ) allocate( PLEPP%dist_locations_i(          n_send       ) )
       if( .not.associated( PLEPP%dist_coords_j)    ) allocate( PLEPP%dist_coords_j(             n_send*ndime ) )
       if( .not.associated( PLEPP%var_ij)           ) allocate( PLEPP%var_ij(            stride, n_send       ) )
       if( .not.associated( PLEPP%var_ji)           ) allocate( PLEPP%var_ji(            stride, n_recv       ) )
       if( .not.associated( PLEPP%interior_list_j)  ) allocate( PLEPP%interior_list_j(           n_recv       ) )
       if(CPLNG_PROPS%ndime_props > 0_ip ) then
          if( .not.associated( PLEPP%dist_props_j)  ) allocate( PLEPP%dist_props_j( n_send * CPLNG_PROPS%ndime_props  ) )
       end if
       PLEPP%n_ij   = n_send
       PLEPP%n_ji   = n_recv
       PLEPP%stride = stride

    endif

  end subroutine plepp_pdn_dynamic_allocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_dynamic_deallocate( PLEPP )

    implicit none

    type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP

    if( ( PLEPP%commij /= -1 .or. PLEPP%commji /= -1 ) .and. initiated ) then

       call commdom_locator_destroy()
       if( associated( PLEPP%tetra_coords_j )  ) deallocate( PLEPP%tetra_coords_j   )
       if( associated( PLEPP%dist_locations_i) ) deallocate( PLEPP%dist_locations_i )
       if( associated( PLEPP%dist_coords_j)    ) deallocate( PLEPP%dist_coords_j    )
       if( associated( PLEPP%var_ij)           ) deallocate( PLEPP%var_ij           )
       if( associated( PLEPP%var_ji)           ) deallocate( PLEPP%var_ji           )
       if( associated( PLEPP%interior_list_j)  ) deallocate( PLEPP%interior_list_j  )
       if(CPLNG_PROPS%ndime_props > 0_ip ) then
          if( associated( PLEPP%dist_props_j)  ) deallocate( PLEPP%dist_props_j     )
       end if
       PLEPP%n_ij   = 0
       PLEPP%n_ji   = 0
       initiated = .false.

    endif

  end subroutine plepp_pdn_dynamic_deallocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_dynamic_get_coupling_properties( PLEPP )

    implicit none

    type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
    integer(ip)                                 :: n_recv, n_send

    if( PLEPP%commij /= -1 .or.  PLEPP%commji /= -1 ) then

       n_recv = 0_ip
       n_send = 0_ip
       call commdom_locator_get_n_dist_points( n_send )
       call commdom_locator_get_n_interior(    n_recv )

       if( INOTMASTER ) then
          PLEPP%dist_locations_i = -1
          PLEPP%dist_coords_j    = -1
          call commdom_locator_get_dist_locations( PLEPP%dist_locations_i )
          call commdom_locator_get_dist_coords(    PLEPP%dist_coords_j    )
          call plepp_pdn_locator_send_nodal_var00( PLEPP%dist_locations_i, &
               PLEPP%dist_coords_j,    &
               n_send, &
               PLEPP%tetra_coords_j )
          if( CPLNG_PROPS%ndime_props > 0_ip ) then
             call commdom_locator_get_dist_props(    PLEPP%dist_props_j    )
          end if
          
          PLEPP%interior_list_j(1:n_recv) = -1
          call commdom_locator_get_interior_list( PLEPP%interior_list_j(1:n_recv) )
          
       endif

    endif

  end subroutine plepp_pdn_dynamic_get_coupling_properties
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_dynamic_reduce_sum(prop_in, prop_out)
    
    implicit none
    
    real(rp),              intent(in ) :: prop_in
    real(rp),              intent(out) :: prop_out
    
    prop_out = huge(1.0_rp)
    if( .not.(PLEPP_CPLNG%commij==-1) ) then
       call commdom_reduce_sum_real(prop_in, prop_out, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commij)
    end if
    if( .not.(PLEPP_CPLNG%commji==-1) ) then
       call commdom_reduce_sum_real(prop_in, prop_out, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commji)
    end if

  end subroutine plepp_pdn_dynamic_reduce_sum

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_dynamic_exchange02(prop_i, prop_j, stride)

    implicit none

    real(rp),    intent(in)  :: prop_i(:,:)
    real(rp),    intent(out) :: prop_j(:,:)
    integer(ip), intent(in)  :: stride

    integer(ip)              :: ii

    if( PLEPP_CPLNG%commij /= -1 .or. PLEPP_CPLNG%commji /= -1 ) then

       if( stride < 1_ip .or. stride > PLEPP_CPLNG%stride ) then
          print *, "plepp_pdn_dynamic_exchange02: stride < PLEPP%stride:", stride, "<",PLEPP_CPLNG%stride
          call runend('EXIT!')
       endif

       PLEPP_CPLNG%var_ij(1_ip:stride,1_ip:PLEPP_CPLNG%n_ij) = 0.0_rp
       if( INOTMASTER ) then
          do ii = 1,stride
             call plepp_pdn_locator_send_nodal_var01(npoin,prop_i(ii,1:npoin),                        & !<---
                  PLEPP_CPLNG%dist_locations_i(  1_ip:PLEPP_CPLNG%n_ij           ), PLEPP_CPLNG%n_ij, &
                  PLEPP_CPLNG%tetra_coords_j(    1_ip:PLEPP_CPLNG%n_ij,1_ip:mnode),                   & !--->
                  PLEPP_CPLNG%var_ij(         ii,1_ip:PLEPP_CPLNG%n_ij           ) )
          enddo
       endif

       PLEPP_CPLNG%var_ji(1_ip:stride,1_ip:PLEPP_CPLNG%n_ji) = 0.0_rp
       call commdom_locator_exchange_double_stride(PLEPP_CPLNG%var_ij, PLEPP_CPLNG%var_ji, stride) !> send, recv
       if( INOTMASTER ) then
          do ii = 1,stride
             prop_j(ii,PLEPP_CPLNG%interior_list_j) = PLEPP_CPLNG%var_ji(ii,1:PLEPP_CPLNG%n_ji)
          enddo
       endif

    endif

  end subroutine plepp_pdn_dynamic_exchange02

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_driver_init(CPLNG)

    use def_master,       only : ID_SOLIDZ
    use def_master,       only : current_code
    use def_coupli,       only : mcoup, RESIDUAL
    use def_kermod,       only : contact_tol
    
    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG

    CPLNG%code_i       =  1_ip           !< CODEi
    CPLNG%module_i     =  ID_SOLIDZ      !< MODULEi
    CPLNG%what_i       =  RESIDUAL       !< 'physical'  or 'numerical' coupling
    
    CPLNG%code_j       =  2_ip           !< CODEj
    CPLNG%module_j     =  ID_SOLIDZ      !< MODULEj
    CPLNG%what_j       = -RESIDUAL       !<---- never
    
    CPLNG%n_dof        =  ndime          !< D.O.F.
     
    PLEPP_CPLNG%tol    =  contact_tol    !< Contact tolerance
    CPLNG_PROPS%ndime_props = ndime      !< props @ vertex_coord_j
    
    CPLNG%current_code = current_code    !< *.dat CODE: ID_CODE
    current_code       = 1_ip            !< trick!! -> PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES
    mcoup              = 0_ip            !< Avoid call cou_turnon!!
    
    !if( IMASTER ) print*, "[commdom_driver_init_contact]"

  end subroutine plepp_pdn_driver_init

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   ???
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_driver_sendrecv(CPLNG, current_when, current_task)

    use def_master,           only: iblok, ittim
    use def_master,           only: modul
    use def_master,           only: nblok
    use def_master,           only: mmodu
    use def_master,           only: ITASK_INIUNK, ITASK_TURNOF
    use def_master,           only: ITASK_TIMSTE, ITASK_ENDSTE
    use def_master,           only: ITASK_BEGZON, ITASK_ENDZON
    use def_master,           only: ITASK_AFTER,  ITASK_BEFORE
    use def_master,           only: ITASK_BEGSTE, ITASK_CONBLK
    use def_master,           only: ID_SOLIDZ

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG
    integer(ip),            intent(in)    :: current_when
    integer(ip),            intent(in)    :: current_task

    integer(ip)                           :: itime
    integer(ip)                           :: now(8_ip)
    integer(ip)                           :: the_time(8_ip,n_times) = -1_ip
    logical(lg)                           :: sendrecv(n_times) = .false.
    !
    !-----------------------------------------------------------------------||---!
    !
    !  ittim: Current time step
    !  itcou: Current global iteration
    !  iblok: Current block
    !  is_when = itcou==micou(iblok)
    !
    now = (/ CPLNG%current_code, iblok, modul, current_task, current_when, ittim, -1_ip, -1_ip /)

    if( nblok>1 ) then
       print *, "[plepp_driver_sendrecv] ", "nblok==1 !!"
       call runend("EXIT!!")
    endif
    !-----------------------------------------------------------| ITERATIONS |---!
    !   +
    !   |_Alya
    !     |_call Domain()
    !     |_call Turnon()
    !     |_call Iniunk()
    !     |_time: do while
    !       |_call Timste()
    !       |_reset: do
    !         |_call Begste()              (+2-)
    !           |_block: do while
    !             |_coupling: do while
    !               |_call Begzon()        (+4-) (+7-)  AFTER<-never INTO the module: into 'coupli/mod_coupling_driver.f90'
    !
    !               |_modules: do while
    !                 |_call Doiter()
    !                 |_call Concou()      (+7-)
    !               |_call Endzon()        (-5+) (+7-) (-12+)
    !             |_call Conblk()          (-6+)
    !       |_call Endste()
    !   |_call Turnof()
    !
    !-----------------------------------------------------------------------||---!
    if( CPLNG%current_code==CPLNG%code_i ) then
       the_time(:,2_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_BEGSTE, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
       the_time(:,4_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_BEGZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
       the_time(:,5_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
       the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_CONBLK, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
       !
       !the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
       the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /) 
       the_time(:,12_ip) = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    end if
    if( CPLNG%current_code==CPLNG%code_j ) then
       the_time(:,2_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGSTE, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
       the_time(:,4_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
       the_time(:,5_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
       the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_CONBLK, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
       !
       the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /) 
       the_time(:,12_ip) = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    end if
    !
    sendrecv = (/ (all(the_time(:,itime)==now), itime=1,n_times) /)
    !
    CNT_SENDRECV = sendrecv                                       
    !
    !                                                                            !
    !-------------------------------------------------------------| -TIMSTE+ |---!
    if( sendrecv(1_ip) ) then

    end if
    !-------------------------------------------------------------| +BEGSTE- |---!
    if( sendrecv(2_ip) ) then
       call plepp_pdn_driver_begste()
    end if
    !-------------------------------------------------------------| -BEGZON+ |---!
    if( sendrecv(4_ip) ) then
       call plepp_pdn_driver_begzon()
    end if
    !-------------------------------------------------------------| -ENDZON+ |---!
    if( sendrecv(5_ip) ) then
       call plepp_pdn_driver_endzon()
    end if
    !-------------------------------------------------------------| +CONBLK- |---!
    if( sendrecv(6_ip) ) then
       call plepp_pdn_driver_conblk()
    end if
    !-------------------------------------------------------------| EXCHANGE |---!
    if( sendrecv(7) .or. sendrecv(12_ip) ) then
       if(      kfl_conta == 1 ) then         ! Unilateral
          if( sendrecv(7_ip) )                      call plepp_pdn_driver_localize_unilateral(CPLNG)
       else if ( kfl_conta == 2 ) then        ! Bilateral
          if( sendrecv(7_ip) .or. sendrecv(12_ip) ) call plepp_pdn_driver_localize_bilateral(CPLNG)  
       end if
       !< AFTER<-never INTO the module !!
       if(current_when==ITASK_AFTER) then
          if( modul==ID_SOLIDZ ) call Solidz( -1_ip )
       endif
       
    endif

  end subroutine plepp_pdn_driver_sendrecv

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   Localization nodes contact for Unilateral cases
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_driver_localize_unilateral(CPLNG)

    use def_master, only : displ
    use def_master, only : ITER_K

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG
    integer(ip)                           :: ipoin
    
    if( INOTMASTER ) then
       !
       ! Unilateral/Neumann contact partd
       !
       do ipoin = 1,npoin
          coord(1:ndime,ipoin) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,ITER_K)
       end do

    end if
    
    ! Dynamic deallocate
    call plepp_pdn_dynamic_deallocate(PLEPP_CPLNG)
    ! Deform mesh
    call plepp_pdn_dynamic_set_mesh(CPLNG%n_dof)

    if( INOTMASTER ) then
       !
       ! Unilateral/Neumann contact parts
       !
       do ipoin = 1,npoin
          coord(1:ndime,ipoin) = coord(1:ndime,ipoin) - displ(1:ndime,ipoin,ITER_K)
       end do
       
    end if

  end subroutine plepp_pdn_driver_localize_unilateral

  !-----------------------------------------------------------------------
  !>
  !> @author  M. Rivero
  !> @date
  !> @brief   Localization nodes contact and deform mesh
  !> @details
  !>          It worked for unilateral only implicit time integration scheme
  !>          not for explicit cases. It is particularly designed for bilateral
  !>          cases.
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_driver_localize_bilateral(CPLNG)

    use def_master, only : displ
    use def_master, only : ITER_K
    use def_coupli, only : RESIDUAL

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG
    integer(ip)                           :: ipoin

    if( INOTMASTER ) then
       !
       ! Unilateral contact part
       !
       if( CNT_SENDRECV(7_ip) ) then   ! solamente para el unilateral
          ! Localizacion i body (A1) with deformed shape
          if( CNT_CPLNG%current_code==CNT_CPLNG%code_i ) then ! Neumann detection (identer)
             do ipoin = 1,npoin
                coord(1:ndime,ipoin) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,ITER_K)
             end do
          end if

          ! Localization j body (A2) with undeformed shape
          if( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) then ! Dirichlet detection (block)
             do ipoin = 1,npoin
                coord(1:ndime,ipoin) = coord(1:ndime,ipoin)
             end do
          end if
       end if
       !
       ! Neumann contact part
       !
       if( CNT_SENDRECV(12_ip) ) then
          ! Localization for both bodies
          do ipoin = 1,npoin
             coord(1:ndime,ipoin) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,ITER_K)
          end do
       end if
    end if
    !
    ! Relaxation (should not be here!!!)
    !
    !if ( CNT_SENDRECV(7_ip ) ) PLEPP_CPLNG%tol =  1.0e-10_rp !< localization tolerance update (defined in mod_commdom_plepp.f90)
    !if ( CNT_SENDRECV(12_ip) ) PLEPP_CPLNG%tol =  1.0e-3_rp  !< localization tolerance update (defined in mod_commdom_plepp.f90)
    !
    ! Deform mesh
    !
    ! Dynamic deallocate
    call plepp_pdn_dynamic_deallocate(PLEPP_CPLNG)
    ! Deform mesh
    call plepp_pdn_dynamic_set_mesh(CPLNG%n_dof)

    if( INOTMASTER ) then
       !
       ! Unilateral contact part
       !
       if( CNT_SENDRECV(7_ip) ) then
          ! Localizacion i body (A1) with deformed shape
          if( CNT_CPLNG%current_code==CNT_CPLNG%code_i ) then ! Neumann detection (identer)
             do ipoin = 1,npoin
                coord(1:ndime,ipoin) = coord(1:ndime,ipoin) - displ(1:ndime,ipoin,ITER_K)
             end do
          end if

          ! Localization j body (A2) with undeformed shape
          if ( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) then ! Dirichlet detection (block)
             do ipoin = 1,npoin
                coord(1:ndime,ipoin) = coord(1:ndime,ipoin)
             end do
          end if

       end if
       !
       ! Neumann contact part
       !
       if( CNT_SENDRECV(12_ip) ) then
          ! Localization for both bodies
          do ipoin = 1,npoin
             coord(1:ndime,ipoin) = coord(1:ndime,ipoin) - displ(1:ndime,ipoin,ITER_K)
          end do
       end if
       !
       ! Residual
       !
       if( resid ) then
       else
          call plepp_pdn_driver_nodes_to_reaction( CPLNG )
          resid = .true.
       endif
       
    end if

  end subroutine plepp_pdn_driver_localize_bilateral

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-26
  !> @brief   Copied from M.Zavala
  !> @details ???
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_driver_nodes_to_reaction( CPLNG )
    
    use def_parame,         only : ip, rp
    use def_master,         only : inotmaster
    use def_domain,         only : npoin, mcono, mcodb,  kfl_codno
    use def_kintyp_solvers, only : soltyp
    use def_master,         only : momod, modul
    use def_coupli,         only : RESIDUAL
    
    implicit none
    
    type(COMMDOM_COUPLING), intent(inout) :: CPLNG
    logical(lg)                           :: sendrecv(5_ip)
    integer(ip)                           :: ipoin, ndofn, ncono
    logical(lg),  pointer                 :: touched(:) => null()
    logical(lg)                           :: codno(mcono)
    logical(lg),  pointer                 :: fixno(:)
    type(soltyp), pointer                 :: solve_sol(:)

    solve_sol => momod(modul) % solve(1_ip:)
    ndofn = solve_sol(1)%ndofn

    if( INOTMASTER ) then
       allocate( touched(npoin) )
       allocate(   fixno(ndofn) )
    else
       allocate( touched(1) )
       allocate(   fixno(1) )
    endif
    
    touched = .false.
    
    if( INOTMASTER ) then
       do ipoin = 1,npoin
          codno(1:mcono) = kfl_codno(1:mcono,ipoin) /= mcodb+1 !< Is it at 'boundary'?
          ncono          = count( codno(1:mcono) , KIND=ip)    !< How many 'codes' are?
          !
          if( ncono>0 ) then  !< Is a                   'boundary'?
             !fixno(1:ndofn) = abs( solve_sol(1)%kfl_fixno(1:ndofn,ipoin) ) == 0    !< Is it fixed? !<2019June17 <GGU>
             if( all(fixno(1:ndofn)) ) touched(ipoin) = .true.
          endif
          !
       enddo
    endif

    CPLNG%current_what = .false.
    CPLNG%current_what(1_ip) = (CPLNG%what_i==RESIDUAL)
    CPLNG%current_what(2_ip) = (CPLNG%what_j==RESIDUAL)
    CPLNG%current_what(3_ip) = CPLNG%current_what(1_ip).and.CPLNG%current_what(2_ip)
    CPLNG%current_what(4_ip) = CPLNG%current_what(1_ip).or. CPLNG%current_what(2_ip)

    sendrecv = .false.
    sendrecv(1_ip) = (CPLNG%code_i==CPLNG%current_code).and.(CPLNG%module_i==modul).and.(CPLNG%current_what(4_ip)).and.INOTMASTER
    sendrecv(2_ip) = (CPLNG%code_j==CPLNG%current_code).and.(CPLNG%module_j==modul).and.(CPLNG%current_what(4_ip)).and.INOTMASTER
    sendrecv(3_ip) = sendrecv(1_ip).and.sendrecv(2_ip)
    sendrecv(4_ip) = sendrecv(1_ip).or. sendrecv(2_ip)
    
    if( sendrecv(4_ip) ) then                    ! OR

       if( sendrecv(2_ip) ) then                 ! CODE_J
          if( CPLNG%current_what(2_ip) ) then    ! BVNAT_J

             call plepp_pdn_allocate_react_bvnat(  momod(modul), modul, .false., .true.)

          else                                    ! REACT_J

             call plepp_pdn_allocate_react_bvnat(  momod(modul), modul, .true., .false.)
             
             call plepp_pdn_set_source_nodes( solve_sol(1)%lpoin_reaction(1:npoin)) ! Mark the nodes where reaction is required
             !solve_sol(1)%lpoin_reaction(1:npoin) = touched(1:npoin)               ! I do not know why, but this doesnt work!! 
             !
             !
             ! LPOIN creation: Sets to true all the lpoin_reaction nodes wihich does not have fixno > 0
             !
             do ipoin = 1,npoin
                if( .not. maxval(solve_sol(1) % kfl_fixno(:,ipoin)) > 0_ip ) solve_sol(1) % lpoin_reaction(1:npoin) = .true.
             end do
             call allocate_block_system( momod(modul), solve_sol(1)%lpoin_reaction(1:npoin) )

          end if

          if(INOTSLAVE) print *, "[commdom_plepp_inivar]", " 'RESIDUALj'", count( solve_sol(1)%lpoin_reaction(1:npoin) , KIND=ip )

       else

          if( sendrecv(1_ip) ) then               !  CODE_I
             if( CPLNG%current_what(1_ip) ) then  ! BVNAT_I

                call plepp_pdn_allocate_react_bvnat(  momod(modul), modul, .false., .true.)
             else                                   ! REACT_I
                
                call plepp_pdn_allocate_react_bvnat(  momod(modul), modul, .true., .false.)
                
                call plepp_pdn_set_source_nodes( solve_sol(1)%lpoin_reaction(1:npoin) ) ! Mark the nodes where reaction is required
                !solve_sol(1)%lpoin_reaction(1:npoin) = touched(1:npoin)                ! I do not know why, but this doesnt work!!
                !
                ! LPOIN creation: Sets to true all the lpoin_reaction nodes wihich does not have fixno > 0
                !
                do ipoin = 1,npoin
                   if( .not. maxval(solve_sol(1) % kfl_fixno(:,ipoin)) > 0_ip ) solve_sol(1) % lpoin_reaction(1:npoin) = .true.
                end do
                call allocate_block_system( momod(modul), solve_sol(1)%lpoin_reaction(1:npoin) )

             end if

             if( INOTMASTER ) print*,"[commdom_plepp_inivar]", " 'RESIDUALi'"

          end if

       end if
       
    end if
    !
    ! Deallocate arrays
    !
    deallocate( touched )
    deallocate( fixno   )

  end subroutine plepp_pdn_driver_nodes_to_reaction


  subroutine plepp_pdn_set_source_nodes(prop_out)
    
    implicit none
    
    logical(lg), intent(inout)  :: prop_out(npoin)
    
    if( PLEPP_CPLNG%commij /= -1 .or. PLEPP_CPLNG%commji /= -1 ) then
       if( INOTMASTER .or. ISEQUEN) then
          prop_out(1_ip:npoin) = .false.
          prop_out(PLEPP_CPLNG%interior_list_j) = .true.
       end if

    endif

  end subroutine plepp_pdn_set_source_nodes

  subroutine plepp_pdn_allocate_react_bvnat(mod_module,module_k,react,bvnat)

    use def_kintyp,         only : ip, rp, lg
    use def_kintyp_solvers, only : soltyp
    use def_kintyp,         only : tymod
    use def_domain,         only : npoin
    use def_postpr,         only : mem_modul
    use mod_memory,         only : memory_alloca
    
    implicit none
    
    integer(ip),  intent(in)    :: module_k
    logical(lg),  intent(in)    :: react,bvnat
    type(tymod),  intent(inout) :: mod_module
    type(soltyp), pointer       :: solve_sol(:)
    integer(ip)                 :: ndofn
    integer(ip)                 :: num_blocks
    integer(ip)                 :: ipoin,i_block,ndofn_block
    
    solve_sol  => mod_module % solve(1:)
    
    if( solve_sol(1) % block_num == 1 .and. solve_sol(1) % kfl_algso /= -999 ) then
       ndofn      = solve_sol(1) % ndofn
       num_blocks = solve_sol(1) % num_blocks
       if( react ) then
          solve_sol(1) % kfl_react = 1
          call memory_alloca(  mem_modul(1:2,module_k), 'SOLVE % REACTION','inivar',solve_sol(1) % reaction, ndofn, npoin)
          do i_block = 2,num_blocks
             ndofn_block = solve_sol(1) % block_dimensions(i_block)
             call memory_alloca(mem_modul(1:2,module_k),'SOLVE % REACTION','inivar',solve_sol(i_block) % reaction,ndofn_block,npoin)
          end do
       end if
       !
       !call memory_alloca(    mem_modul(1:2,module_k), 'SOLVE_SOL(1) % LPOIN_REACTION','inivar',solve_sol(1) % lpoin_reaction,npoin)
       !solve_sol(1) % lpoin_reaction(1:npoin) = .false.
       !
       if( bvnat ) then
          solve_sol(1) % kfl_bvnat = 1
          call memory_alloca(  mem_modul(1:2,module_k), 'SOLVE % BVNAT','inivar',solve_sol(1) % block_array(1) % bvnat,ndofn, npoin)
          solve_sol(1) % bvnat => solve_sol(1) % block_array(1) % bvnat
          do i_block = 2,num_blocks
             ndofn_block = solve_sol(1) % block_dimensions(i_block)
             call memory_alloca(mem_modul(1:2,module_k),'SOLVE % BVNAT','inivar',solve_sol(1) % block_array(i_block) % bvnat, ndofn_block, npoin)
          end do
       end if
       
       call memory_alloca(mem_modul(1:2,module_k), 'SOLVE_SOL % LPOIN_REACTION','inivar', solve_sol(1) % lpoin_reaction, npoin)
       do ipoin = 1,npoin
          solve_sol(1) % lpoin_reaction(ipoin) = .false.
       end do
       
    end if
    
  end subroutine plepp_pdn_allocate_react_bvnat

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------
  
  subroutine plepp_pdn_driver_exchange02(CPLNG)

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG

    call plepp_pdn_dynamic_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
         CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
         CPLNG%n_dof)

  end subroutine plepp_pdn_driver_exchange02

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_driver_memall(CPLNG)

    use def_domain, only: npoin

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG

    CPLNG%n_pts = 0_ip
    if( INOTMASTER ) then
       CPLNG%n_pts = npoin
       !CPLNG%n_pts = CPLNG_PROPS%npoin
    end if

    if( .not.associated(CPLNG%var_ij) ) then
       allocate( CPLNG%var_ij(CPLNG%n_dof,CPLNG%n_pts) )
    else
       call runend('[commdom_alya_memall] ERROR: code_i==code_j!!')
    end if
    if( .not.associated(CPLNG%var_ji) ) then
       allocate( CPLNG%var_ji(CPLNG%n_dof,CPLNG%n_pts) )
    end if

    CPLNG%var_ji(1:CPLNG%n_dof,1:CPLNG%n_pts) = 0.0_rp
    CPLNG%var_ij(1:CPLNG%n_dof,1:CPLNG%n_pts) = 0.0_rp

    !if( IMASTER .or. ISEQUEN ) print*, "plepp_pdn_driver_memall"

  end subroutine plepp_pdn_driver_memall

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_driver_begste()

    use def_master,   only : iblok
    use def_master,   only : ITASK_TIMSTE
    use def_master,   only : dtinv, cutim, dtime
    use def_coupli,   only : coupling_driver_iteration
    use def_coupli,   only : coupling_driver_number_couplings
    use def_coupli,   only : coupling_driver_max_iteration
    use def_coupli,   only : max_block_cou
    use mod_messages, only : livinf
    
    implicit none

    coupling_driver_iteration(1:max_block_cou) = 0_ip
    coupling_driver_number_couplings(iblok)    = 1_ip
    coupling_driver_max_iteration(iblok)       = 1_ip

    call plepp_pdn_compare_dtinv(dtinv) 
    cutim  = cutim - dtime
    call setgts(ITASK_TIMSTE)
    call livinf(201_ip, ' ',1_ip)

  end subroutine plepp_pdn_driver_begste

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_driver_begzon()

    use def_master,    only : iblok
    use def_master,    only : mmodu
    use def_master,    only : lmord
    use def_master,    only : itinn
    use def_coupli,    only : coupling_driver_iteration
    use def_coupli,    only : coupling_driver_number_couplings
    use mod_messages,  only : livinf

    implicit none

    integer(ip) :: iorde,imodu

    if( coupling_driver_number_couplings(iblok) /= 0 .and. &
        coupling_driver_iteration(iblok) == 0 ) then
       call livinf(-6_ip,'ZONAL COUPLING FOR BLOCK ', iblok)
    end if
    !
    ! Put inner iterations to zero
    !
    do iorde = 1,mmodu
       imodu = lmord(iorde,iblok)
       itinn(imodu) = 0
    end do
    !
    ! Iteration counter
    !
    coupling_driver_iteration(iblok) = coupling_driver_iteration(iblok) + 1

  end subroutine plepp_pdn_driver_begzon

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_driver_endzon()

    use def_master,  only : iblok
    use def_master,  only : kfl_gocou
    use def_coupli,  only : coupling_driver_iteration
    use def_coupli,  only : coupling_driver_number_couplings
    use def_coupli,  only : coupling_driver_max_iteration
    use def_coupli,  only : kfl_gozon

    implicit none

    if( coupling_driver_number_couplings(iblok) /= 0 ) then
       
       kfl_gozon = 0
       if( coupling_driver_iteration(iblok) >= coupling_driver_max_iteration(iblok) ) then
          kfl_gozon = 0 
       else
          kfl_gozon = 1
       endif
       if( kfl_gozon == 1 ) kfl_gocou = 1
    else
       
       kfl_gozon = 0
       
    end if

  end subroutine plepp_pdn_driver_endzon

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam
  !> @date    2022-05-25
  !> @brief   
  !> @details Copied from JM Zavala-Ake
  !> 
  !-----------------------------------------------------------------------

  subroutine plepp_pdn_driver_conblk()

    use def_master,   only : iblok
    use def_coupli,   only : coupling_driver_number_couplings
    use def_coupli,   only : coupling_driver_iteration
    use mod_messages, only : livinf

    implicit none
    !
    ! Write message if block is in a coupling loop
    !
    if( coupling_driver_number_couplings(iblok) /= 0 ) then
       call livinf(-13_ip,'END ZONAL COUPLING: ', coupling_driver_iteration(iblok) )
    end if

  end subroutine plepp_pdn_driver_conblk

#endif

end module mod_plepp_pdn_contact


