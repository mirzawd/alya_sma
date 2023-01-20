!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    ker_extrap_boundary.f90
!> @author  Herbert Owen
!> @brief    
!> @details - Initialization to obtain value (for example shear stress) from closest point at the boundary
!!          - needed for mixing length turbulence moder for non equilibrium 2 layer wall model
!!          - generates: ptb_to_use, wallcoupling_extr_boun and is_interior
!!          - First for each 'interior' node get lninv_loc of the nearest boundary node. Also get the coordinates of the boundary point. 
!!          - 'interior' here stands for not belonging to the boundaries of interest wich can be selectected with RESTR in the ker.dat file.
!!          - Sort according to lninv_loc of the boundary node. this allows to groups all interior points that belong to the same boundary node.
!!          - This information is saved in ptb_to_use.
!!          - Second part obtains wallcoupling_extr_boun that will allow to obtain whatever value is needed(for example velocity) at the boundary points. 
!> @} 
!-----------------------------------------------------------------------

subroutine ker_extrap_boundary()
  use def_domain
  use def_kintyp,                     only     : ip,rp,lg
  use def_master,                     only     : kfl_paral, current_code,INOTMASTER, zeror, ID_TEMPER,ID_NASTIN
  use def_master,                     only     : kfl_modul, lzone , modul, mem_modul, lninv_loc
  use def_kermod,                     only     : wallcoupling_extr_boun, ptb_to_use, is_interior, nrestr_2_codno, krestr_2_codno, kfl_twola_ker
  use mod_interpolation,              only     : COU_GET_INTERPOLATE_POINTS_VALUES

  use mod_couplings,                  only     : COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_couplings_communications,   only     : COU_GENERATE_LOCAL_TRANSMISSION_MATRICES
  use mod_couplings_communications,   only     : COU_PARALLELIZE_TRANSMISSION_MATRICES
  use mod_coupling_memory,            only     : cou_initialization  
  use def_coupli,                     only     : NEAREST_BOUNDARY_NODE
  use def_coupli,                     only     : typ_color_coupling
  use mod_parall,                     only     : par_code_zone_subd_to_color,PAR_MY_CODE_RANK
  use mod_parall,                     only     : PAR_COMM_COLOR
  use mod_memory,                     only     : memory_alloca,memory_deallo
  use mod_maths,                      only     : maths_heap_sort
  use def_domain,                     only     : mcono,kfl_codno
  use mod_communications,             only     : PAR_BARRIER


  implicit none

  integer(ip)              :: ipoin,ii,ierr,icono
  integer(ip)              :: icolo,jcolo,izone,size1,size2

  type(typ_color_coupling) :: wallcoupling_lninv_loc  ! used to extrapolate shear stress from boundaries for two layer model

  ! see what to do if master -- see from where we call this

  integer(ip)               :: iauxi,np_int,kount,kurrent
  integer(ip) , pointer     :: lninv_loc_boun(:)
  real(rp)    , pointer     :: dlninv_loc(:)      ! to save the fact that NO specific subroutine for call to COU_GET_INTERPOLATE_POINTS_VALUES
  real(rp)    , pointer     :: dlninv_loc_boun(:)
  integer(ip) , pointer     :: aux_copy(:)
  integer(ip) , pointer     :: index_table(:)
  real(rp)    , pointer     :: coord_int(:,:)
  real(rp)    , pointer     :: coord_boun(:,:)
  real(rp)    , pointer     :: coord_ptb(:,:)
  logical(lg) , pointer     :: selected_boundary_node(:)

  nullify(is_interior)
  nullify(lninv_loc_boun)
  nullify(dlninv_loc)
  nullify(dlninv_loc_boun)
  nullify(aux_copy)
  nullify(index_table)
  nullify(coord_int)
  nullify(coord_boun)
  nullify(coord_ptb)
  nullify(ptb_to_use)
  nullify(selected_boundary_node)

  kount = -1_ip

  call PAR_BARRIER()


!  if ( kfl_twola_ker > 0 ) then 
  if ( kfl_twola_ker == 2 ) then 
     !
     ! Mark 'interior' - non wall law boundaries -- perhaps we could also want wall -- See what georgios is using
     !
     if( INOTMASTER ) then   ! OJO aca creo que la estoy cagando

        call memory_alloca(mem_modul(1:2,modul),'is_interior',           'ker_extrap_boundary',is_interior            ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'SELECTED_BOUNDARY_NODE','ker_extrap_boundary',selected_boundary_node ,npoin)
        !
        ! Restrict source nodes only to some codes
        ! Actually I could work only with selected_boundary_node and avoid is_interior but i would need to reorder a bit
        !
        is_interior = 1_ip
        selected_boundary_node = .false.
        do ipoin = 1,npoin
           do ii = 1,nrestr_2_codno
              do icono = 1, mcono
                 if ( kfl_codno(icono,ipoin) == krestr_2_codno(ii) ) then
                    selected_boundary_node(ipoin) = .true.
                    is_interior(ipoin) = 0_ip
                 end if
              end do
           end do
        end do
        !
        ! Count
        !
        np_int = 0_ip
        do ipoin=1,npoin
           if ( is_interior(ipoin) == 1_ip ) np_int = np_int + 1_ip
        end do
        call memory_alloca(mem_modul(1:2,modul),'COORD_INT',      'ker_extrap_boundary',coord_int,      ndime,np_int)
        call memory_alloca(mem_modul(1:2,modul),'coord_boun',     'ker_extrap_boundary',coord_boun,     ndime,np_int)
        call memory_alloca(mem_modul(1:2,modul),'lninv_loc_boun', 'ker_extrap_boundary',lninv_loc_boun, np_int)
        call memory_alloca(mem_modul(1:2,modul),'dlninv_loc'     ,'ker_extrap_boundary',dlninv_loc     ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'dlninv_loc_boun','ker_extrap_boundary',dlninv_loc_boun,np_int)
        call memory_alloca(mem_modul(1:2,modul),'aux_copy',       'ker_extrap_boundary',aux_copy,       np_int)
        call memory_alloca(mem_modul(1:2,modul),'index_table',    'ker_extrap_boundary',index_table,    np_int)
        !
        ! First for all interior nodes find the lninv_loc of the closest boundary node  - I call them interior but it is actually all non wall law boundary nodes
        !
        !  put coordinates of all interior nodes in a coord_int(ndime,np_int)
        np_int = 0_ip
        do ipoin=1,npoin
           if ( is_interior(ipoin) == 1_ip ) then
              np_int = np_int + 1_ip
              coord_int(:,np_int) = coord(:,ipoin)
           end if
        end do

     else   ! MASTER
        np_int = 0_ip   ! temporary solution else line wallcoupling_lninv_loc % wet % npoin_wet    = np_int   got complaint in debug
     end if ! Not MASTER
     !
     ! Preliminaries wallcoupling_lninv_loc
     !
     !
     ! Intialize interpolation
     !
     icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip) 
     jcolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip) 

     !
     ! Basic information needed to fill the coupling structures
     !
     call cou_initialization(wallcoupling_lninv_loc)

     wallcoupling_lninv_loc % number             = 1020_ip               ! Coupling number
     wallcoupling_lninv_loc % itype = NEAREST_BOUNDARY_NODE
     wallcoupling_lninv_loc % kfl_toda_costa     = 0_ip       
     wallcoupling_lninv_loc % color_source       = jcolo                 ! source color
     wallcoupling_lninv_loc % color_target       = icolo                 ! target color   
     wallcoupling_lninv_loc % zone_source        = 0_ip                 ! current_zone
     wallcoupling_lninv_loc % zone_target        = 0_ip                 ! current_zone
     wallcoupling_lninv_loc % wet % npoin_wet    = np_int
     wallcoupling_lninv_loc % wet % number_wet_points  = 0_ip

     if( associated(coord_int) ) then
        size1 = size(coord_int,1)
        size2 = size(coord_int,2)
        allocate(wallcoupling_lninv_loc % wet % coord_wet(size1,size2))
        wallcoupling_lninv_loc % wet % coord_wet = coord_int
     end if
     if(kfl_twola_ker==1) print*,'4ker_extrap_boundary:',kfl_paral

     call COU_INIT_INTERPOLATE_POINTS_VALUES(coord_int,icolo,jcolo,wallcoupling_lninv_loc,CANDIDATE_SOURCE_NODES=selected_boundary_node)
     if(kfl_twola_ker==1) print*,'5ker_extrap_boundary:',kfl_paral

     !
     ! check if everything was right
     !
     ierr = 0
     do ii = 1,kount
        if( wallcoupling_lninv_loc % geome % status(ii) == 0 ) then
           print*,ii,PAR_MY_CODE_RANK,' is lost with coordinates',coord_int(:,ii),'kfl_paral',kfl_paral,'wallcoupling_lninv_loc'
           ierr = ierr + 1
        end if
     end do
     if ( ierr/=0 ) call runend('QUE PASA AQUI wallcoupling_lninv_loc ')
     !
     ! Communication matrix creation
     !
     call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(wallcoupling_lninv_loc)
     !
     ! get lninv_loc_boun  & coord_boun
     ! ojo con lninv_loc me da error #6285: There is no matching specific subroutine for this generic subroutine call.   [COU_GET_INTERPOLATE_POINTS_VALUES]
     ! call COU_GET_INTERPOLATE_POINTS_VALUES(lninv_loc, lninv_loc_boun, wallcoupling_lninv_loc)
     !
     dlninv_loc = dble(lninv_loc(1:npoin))
     dlninv_loc_boun = 0.0_rp
     call COU_GET_INTERPOLATE_POINTS_VALUES(dlninv_loc, dlninv_loc_boun, wallcoupling_lninv_loc)    ! lninv_loc_boun(np_int)    -- global numbering of the closest boundary point

     lninv_loc_boun = nint(dlninv_loc_boun)
     call COU_GET_INTERPOLATE_POINTS_VALUES(coord,      coord_boun,      wallcoupling_lninv_loc)    ! coord_boun(np_int)    -- coordinates of the closest boundary point
     if( INOTMASTER ) then   

        !
        ! sort lninv_loc_boun
        !
        do ipoin=1,np_int     ! Order by lninv_loc_boun
           index_table(ipoin) = ipoin
           aux_copy(ipoin) = lninv_loc_boun(ipoin)
        end do
        iauxi = np_int   ! because it is intent inout in maths_heap_sort
        call maths_heap_sort(2_ip,iauxi,aux_copy,' ',index_table)  ! 2 increasing

        !
        !  Obtain kount - number of boundary point that are needed to extrapolate
        !
        kount = 1
        kurrent = aux_copy(1)
        do ipoin=2,np_int
           if (aux_copy(ipoin) /= kurrent) then
              kount = kount + 1
              kurrent = aux_copy(ipoin)
           end if
        end do

        call memory_alloca(mem_modul(1:2,modul),'coord_ptb',     'ker_extrap_boundary',coord_ptb,     ndime,kount)

        !
        ! Obtain ptb_to_use that gives for each interior node the boundary node from which it gets its value
        !        coord_ptb -- the coordinates of the boundary point from which values are extrapolated
        !
        kount = 1
        ipoin = 1
        call memory_alloca(mem_modul(1:2,modul),'ptb_to_use',   'ker_extrap_boundary',ptb_to_use,   np_int)
        ptb_to_use(index_table(ipoin)) = kount
        coord_ptb(:,kount) = coord_boun(:,index_table(ipoin))
        kurrent = aux_copy(1)
        do ipoin=2,np_int
           if (aux_copy(1) /= kurrent) then
              kount = kount + 1
              kurrent = aux_copy(ipoin)
              coord_ptb(:,kount) = coord_boun(:,index_table(ipoin))
           end if
           ptb_to_use(index_table(ipoin)) = kount
        end do

        !
        ! Testing 
        !
        do ipoin=1,np_int
           write(*,'(i9,15(1x,e10.3))') &
                'lninv_loc(index_table(ipoin)),coord(:,index_table(ipoin)),coord_ptb(:,ptb_to_use(index_table(ipoin)) )', &
                lninv_loc(index_table(ipoin)),coord(:,index_table(ipoin)),coord_ptb(:,ptb_to_use(index_table(ipoin)) )
        end do
        !   
        ! wallcoupling_lninv_loc can now be eliminated    
        !

        !  SECOND PART -------------------------------------------------------------------------------------------------------------------------

        !
        ! prliminaries wallcoupling_extr_boun
        !
        !
        ! Intialize interpolation
        !
     end if ! NOTMASTER

     if (1==2) then  !para probar solo al primera paret por ahora

        if (kfl_modul(ID_NASTIN) == 1) then
           izone = lzone(ID_NASTIN)
        else if (kfl_modul(ID_TEMPER) == 1) then
           izone =lzone(ID_TEMPER)
        end if
        icolo = par_code_zone_subd_to_color(current_code,izone,0_ip)     ! icolo and jcolo I do not understan well 
        jcolo = par_code_zone_subd_to_color(current_code,izone,0_ip)     ! I just put them following nsi_velobl
        !
        ! Basic information needed to fill the coupling structures
        !
        call cou_initialization(wallcoupling_extr_boun)
        wallcoupling_extr_boun % number             = 1021_ip               ! Coupling number
        wallcoupling_extr_boun % itype = NEAREST_BOUNDARY_NODE
        wallcoupling_extr_boun % kfl_toda_costa     = 0_ip       
        wallcoupling_extr_boun % color_source       = icolo                 ! source color
        wallcoupling_extr_boun % color_target       = jcolo                 ! target color   
        wallcoupling_extr_boun % zone_source        = izone                 ! current_zone
        wallcoupling_extr_boun % zone_target        = izone                 ! current_zone
        wallcoupling_extr_boun % wet % npoin_wet     = kount
        wallcoupling_extr_boun % wet % number_wet_points  = 0_ip
        wallcoupling_extr_boun % commd % PAR_COMM_WORLD   = PAR_COMM_COLOR(icolo,jcolo)

        if( associated(coord_ptb) ) then
           size1 = size(coord_ptb,1)
           size2 = size(coord_ptb,2)
           allocate(wallcoupling_extr_boun % wet % coord_wet(size1,size2))
           wallcoupling_extr_boun % wet % coord_wet = coord_ptb
        end if
        call COU_INIT_INTERPOLATE_POINTS_VALUES(coord_ptb,icolo,jcolo,wallcoupling_extr_boun,CANDIDATE_SOURCE_NODES=selected_boundary_node)

        call memory_deallo(mem_modul(1:2,modul),'SELECTED_BOUNDARY_NODE','ker_extrap_boundary',selected_boundary_node)
        !
        ! check if everything was right
        !
        ierr = 0
        do ii = 1,kount
           if( wallcoupling_extr_boun % geome % status(ii) == 0 ) then
              print*,ii,PAR_MY_CODE_RANK,' is lost with coordinates',coord_ptb(:,ii),'kfl_paral',kfl_paral
              ierr = ierr + 1
           end if
        end do
        if ( ierr/=0 ) call runend('QUE PASA AQUI')
        !
        ! Communication matrix creation
        !   
        call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(wallcoupling_extr_boun)

     end if

  end if ! kfl_twola_ker == 1


end subroutine ker_extrap_boundary

