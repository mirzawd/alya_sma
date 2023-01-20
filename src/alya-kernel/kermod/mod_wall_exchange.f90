!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_wall_exchange.f90
!> @author  Herbert Owen  & Matias Avila & houzeaux
!> @date    2022-01-17
!> @brief   Wall exchange subroutines
!> @details Wall exchange subroutines
!-----------------------------------------------------------------------

module mod_wall_exchange

  use def_kintyp,                   only : ip,rp
  use def_master,                   only : kfl_paral, current_code,INOTMASTER,IPARALL
  use def_master,                   only : zeror, ID_TEMPER,ID_NASTIN
  use def_master,                   only : kfl_modul, lzone, IMASTER
  use def_master,                   only : modul,mem_modul,momod
  use def_master,                   only : INOTSLAVE
  use def_master,                   only : lun_tempo
  use def_domain,                   only : nnode,ltype
  use def_domain,                   only : ndime,ngaus
  use def_domain,                   only : lnods,lnodb
  use def_domain,                   only : ltype,ltypb
  use def_domain,                   only : lelbo,elmar
  use def_domain,                   only : kfl_codbo
  use def_domain,                   only : ywalb,ndimb
  use def_domain,                   only : mnode,mnodb
  use def_domain,                   only : memor_dom
  use def_domain,                   only : nboun,coord
  use def_domain,                   only : meshe,mgaub
  use def_kermod,                   only : kfl_delta, dexlo_ker, velel_ker
  use def_kermod,                   only : temel_ker, lexlo_ker, wallcoupling_waexl
  use def_kermod,                   only : kfl_waexl_ker,kfl_waexl_imp_ker
  use def_kermod,                   only : kfl_boexc_ker,interp_waexlo,shape_waexl
  use def_kermod,                   only : search_waexlo_seq
  use def_kermod,                   only : search_waexlo_par
  use mod_elsest,                   only : elsest_host_element
  use def_kermod,                   only : delta_dom,ielse,relse,ndivi
  use mod_messages,                 only : messages_live
  use mod_communications,           only : PAR_BARRIER
  use mod_communications,           only : PAR_GATHER
  use mod_communications,           only : PAR_GATHERV
  use mod_communications,           only : PAR_BARRIER
  use mod_communications,           only : PAR_MIN
  use mod_memory,                   only : memory_alloca
  use mod_memory,                   only : memory_deallo
  use mod_memory,                   only : memory_size
  use def_interpolation_method,     only : interpolation
  use def_interpolation_method,     only : INT_ELEMENT_INTERPOLATION
  use mod_parall,                   only : PAR_COMM_MY_CODE
  use mod_parall,                   only : PAR_CODE_SIZE
  use def_master,                   only : npart
  use mod_communications_global,    only : PAR_SUM
  use mod_par_additional_arrays,    only : par_bounding_box
  use mod_bouder,                   only : bouder
  use mod_iofile,                   only : iofile
  use mod_iofile,                   only : iofile_flush_unit
  use mod_par_output_partition,     only : par_output_search_timings
  use mod_strings,                  only : integer_to_string
  use mod_strings,                  only : real_to_string
  implicit none

  private

  interface ker_waexlo_getval
     module procedure &
          ker_waexlo_getval_22,&
          ker_waexlo_getval_21,&
          ker_waexlo_getval_32
  end interface ker_waexlo_getval

  interface ker_waexlo_getder
     module procedure &
          ker_waexlo_getder_22
  end interface ker_waexlo_getder

  public :: ker_waexlo
  public :: ker_adapel
  public :: ker_waexlo_getval
  public :: ker_waexlo_getder

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  Herbert Owen  & Matias Avila & houzeaux
  !> @date    2022-01-17
  !> @brief   Preliminary operations to apply the wall law at the 'Exchange Location' 
  !> @details - See Bodart and Larsson - Wall modeled large eddy simulation in complex ... 2011
  !>          - Obtain the exchange location for each wall law boundary gauss point
  !>          - Initialize interpolation - PAR_INIT_INTERPOLATE_POINTS_VALUES
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_waexlo()

    integer(ip)              :: ielem,inode,ipoin,kount,ii,ierr
    integer(ip)              :: pnode,pgaus,iboun,igaub,inodb
    integer(ip)              :: pelty,pblty,pnodb,pgaub
    real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode),gbcod(ndime)
    real(rp)                 :: baloc(ndime,ndime),eucta
    real(rp)                 :: shapf(64),deriv(64*3),coloc(3), dista, dist_aux
    character(200)           :: aux_char
    real(rp)                 :: time1,time2
    integer(ip),   pointer   :: ielem_waexl(:)   
    real(rp),      pointer   :: cooel_ker(:,:)  
    real(rp),      pointer   :: bobox(:,:,:)
    real(rp),      pointer   :: subox(:,:,:)

    nullify(bobox)
    nullify(subox)
    nullify(ielem_waexl)
    nullify(cooel_ker)

    if( kfl_waexl_ker /= 0_ip ) then 

       call cputim(time1) 
       call messages_live('WALL EXCHANGE','START SECTION')

       !-------------------------------------------------------------------
       !
       ! Define search strategies
       !
       !-------------------------------------------------------------------
       !
       ! Bounding boxes
       !
       call messages_live('DEFINE SEARCH STRATEGY')
       call meshe(ndivi) % element_bb  (bobox,MEMORY_COUNTER=memor_dom)
       call par_bounding_box           (PAR_COMM_MY_CODE,bobox,subox,MEMORY_COUNTER=memor_dom)
       !
       ! Set(fill) search method 
       !
       call search_waexlo_seq % set    (BOBOX=bobox,NAME='WAEXLO')
       call search_waexlo_par % set    (BOBOX=subox,NAME='WAEXLO')

       !print*,'A=',kfl_paral,search_waexlo_seq % type,search_waexlo_par % type

       !-------------------------------------------------------------------
       !
       ! To guarantee that all exchange location points will be found.
       ! ywalb(iboun) = ywalb(iboun) * fact(itest) 
       !
       !-------------------------------------------------------------------

       if( kfl_delta == 1 ) call ker_adapel()

       !-------------------------------------------------------------------
       !
       ! Count elegible boundaries KOUNT and allocate memory
       !
       !-------------------------------------------------------------------

       kount = 0_ip
       if( delta_dom > zeror .or. kfl_delta == 1 ) then
          do iboun = 1,nboun
             if( kfl_boexc_ker(kfl_codbo(iboun)) == 1 ) then
                pblty = ltypb(iboun) 
                kount = kount + ngaus(pblty)
             end if
          end do
       end if
       !
       ! Permanent memory
       !
       call memory_alloca(mem_modul(1:2,modul),'VELEL_KER'  ,'ker_waexlo',velel_ker  ,ndime,kount)
       call memory_alloca(mem_modul(1:2,modul),'LEXLO_KER'  ,'ker_waexlo',lexlo_ker  ,mgaub,nboun)
       call memory_alloca(mem_modul(1:2,modul),'SHAPE_WAEXL','ker_waexlo',shape_waexl,mnode,kount)
       if (kfl_modul(ID_TEMPER) == 1) call memory_alloca(mem_modul(1:2,modul),'TEMEL_KER','ker_waexlo',temel_ker,kount) 
       !
       ! Local memory
       !
       call memory_alloca(mem_modul(1:2,modul),'IELEM_WAEXL','ker_waexlo',ielem_waexl,kount)
       call memory_alloca(mem_modul(1:2,modul),'COOEL_KER'  ,'ker_waexlo',cooel_ker,ndime,kount)

       !-------------------------------------------------------------------
       !
       ! Loop over boundaries to compute COOEL_KER 
       !
       !-------------------------------------------------------------------

       kount = 0_ip
       dist_aux = dexlo_ker

       if( delta_dom > zeror .or. kfl_delta == 1 )  then

          boundaries: do iboun = 1,nboun

             if( kfl_boexc_ker(kfl_codbo(iboun)) == 1 ) then
                !
                ! Variable wall distance (use ywalb instead of dexlo_ker)
                !
                if( kfl_delta == 1 ) dist_aux = ywalb(iboun)     
                !
                ! Element properties and dimensions
                !
                pblty = ltypb(iboun) 
                pnodb = nnode(pblty)
                ielem = lelbo(iboun)
                pelty = ltype(ielem)
                pnode = nnode(pelty)
                pgaub = ngaus(pblty) 
                pgaus = ngaus(pelty)
                !
                ! Gather operations: ELCOD, BOCOD
                !
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)             
                end do
                do inodb = 1,pnodb 
                   ipoin = lnodb(inodb,iboun)
                   bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
                end do

                gauss_points: do igaub = 1,pgaub
                   !
                   ! Obtain normal (baloc(:,ndime) to the surface (following nsi_bouset)
                   !                 
                   call bouder(&
                        pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                        bocod,baloc,eucta)                                   ! and Jacobian
                   call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                   kount = kount + 1_ip
                   lexlo_ker(igaub,iboun) = kount

                   gbcod = 0.0_rp
                   do inodb = 1,pnodb
                      gbcod(1:ndime) = gbcod(1:ndime) &
                           + elmar(pblty) % shape(inodb,igaub) * bocod(1:ndime,inodb)
                   end do
                   !
                   ! for the moment I will set dexlo_ker a const
                   ! value. Later we can set some constant times
                   ! de 1st elemente normal size or even more elaborate opts
                   !
                   cooel_ker(1:ndime,kount) = gbcod(1:ndime) - dist_aux * baloc(1:ndime,ndime) 

                end do gauss_points

             end if

          end do boundaries

       end if

       !-------------------------------------------------------------------
       !
       ! Interpolation strategy
       !
       !-------------------------------------------------------------------

       call messages_live('COMPUTE INTERPOLATION')

       call interp_waexlo % init       ()
       call interp_waexlo % input      (search_waexlo_seq % method,search_waexlo_par % method,&
            &                          COMM=PAR_COMM_MY_CODE,&
            &                          INTERPOLATION_METHOD=INT_ELEMENT_INTERPOLATION,&
            &                          NAME='WALL EXCHANGE',&
            &                          DERIVATIVES=.true.)
       call interp_waexlo % preprocess (cooel_ker,meshe(ndivi))
       call interp_waexlo % output     (UNIT=momod(modul) % lun_outpu,HEADER='        ')

       call memory_deallo(memor_dom,'BOBOX','ker_waexlo',bobox)
       call memory_deallo(memor_dom,'SUBOX','ker_waexlo',subox)

       call par_output_search_timings(search_waexlo_seq % method)

       !-------------------------------------------------------------------
       !
       ! Check everything's ok
       !
       !-------------------------------------------------------------------

       ierr = 0
       loop_ii: do ii = 1,kount
          if( .not. interp_waexlo % found(ii) ) then
             ierr = ierr + 1
             exit loop_ii
          end if
       end do loop_ii
       call PAR_SUM(ierr)
       if( ierr /= 0 ) call runend('WALL EXCHANGE: POINTS NOT FOUND')

       !-------------------------------------------------------------------
       !
       ! Implicit case:
       ! Obtain shapes and ielem for each (iboun, igaub)
       !
       !-------------------------------------------------------------------

       if( kfl_waexl_imp_ker /= 0 ) then
          !
          ! elsest does not work in parallel - for the implicit case, since you always search in the first element,
          ! there is no problem
          !
          do ii = 1,kount
             call elsest_host_element(&
                  ielse,relse,1_ip,meshe(ndivi),cooel_ker(:,ii),ielem,&
                  shapf,deriv,coloc,dista)
             if( ielem > 0 ) then
                pelty = abs(ltype(ielem))
                pnode = nnode(pelty)
                ielem_waexl(ii) = ielem
                shape_waexl(1:pnode,ii) = shapf(1:pnode)
             else
                write(aux_char,*)'ker_waexlo:element not found for cooel_ker(:,ii)',cooel_ker(:,ii),ii,kfl_paral
                call runend(aux_char)
             end if
          end do
          ! if implicit exchange location, confirm point is inside first element

          do iboun =1, nboun
             pblty = ltypb(iboun) 
             pnodb = nnode(pblty)
             ielem = lelbo(iboun)
             do igaub = 1,pgaub
                ii =  lexlo_ker(igaub,iboun)
                if (ielem_waexl(ii)/=ielem) then
                   print *, 'ker_waexlo:bad_element, iproce, ielem,connec_elem:', kfl_paral,ielem_waexl(ii),ielem 
                   call runend('ker_waexlo: could not perform exch location implicit because point outside first element')     
                end if
             end do
          end do
       end if

       call memory_deallo(mem_modul(1:2,modul),'IELEM_WAEXL','ker_waexlo',ielem_waexl)
       call memory_deallo(mem_modul(1:2,modul),'COOEL_KER'  ,'ker_waexlo',cooel_ker)

       call cputim(time2) ; time2 = time2 - time1
       call messages_live('TOTAL CPU TIME= '//real_to_string(time2,'(e13.6)'))
       call messages_live('WALL EXCHANGE','END SECTION')

    end if

  end subroutine ker_waexlo

  !-----------------------------------------------------------------------
  !> 
  !> @author  Herbert Owen and houzeaux
  !> @date    2022-01-17
  !> @brief   ADAPt Exchange Location -modifies ywalb(iboun) = ywalb(iboun) * fact(itest) ; to guarantee that all exchange location points will be found.
  !> @details - Adapt the distance used in exchange location so that all points are found
  !>          - This will be necesary for complex geometries moreover it will give robustness in case elsest fails by slightly moving the point.
  !>          - It only works with variable wall distance.
  !>          - The points that are modified can be seen in Paraview. See below.
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_adapel()

    integer(ip)                          :: ielem,inode,ipoin,kount
    integer(ip)                          :: pnode,pgaus,iboun,igaub,inodb
    integer(ip)                          :: pelty,pblty,pnodb,pgaub
    real(rp)                             :: bocod(ndime,mnodb)
    real(rp)                             :: elcod(ndime,mnode),gbcod(ndime)
    real(rp)                             :: baloc(ndime,ndime),eucta
    integer(ip)                          :: kfl_all_found
    real(rp)                             :: dist_aux
    integer(ip),     parameter           :: ntest = 11
    real(rp)                             :: fact(ntest)
    real(rp),        pointer             :: xcoor(:,:)
    integer(ip),     pointer             :: missing_find(:)
    integer(ip),     pointer             :: level_found(:)
    integer(ip)                          :: itest,i,ii,ipart

    type(interpolation)                  :: fake_waexlo                ! Wall exchange iterpolation  
    !
    ! for visualization - Ideas borrowed from outstl
    !
    character(150)                       :: fil_tempo
    real(rp),        pointer             :: aux_vec(:,:)
    real(rp),        pointer             :: aux_vec_gat(:)   ! in outstl xstl_gat(:)  but I belive it is cleaner to have it (:,:) - but it does not work
    integer(ip),     parameter           :: n1_aux_vec = 9_ip
    integer(4)                           :: nauxv4,nauxv4_tot    
    integer(4),      pointer             :: nauxv4_gat(:)
    integer(ip)                          :: kount_not_found_in_first
    real(rp)                             :: time1,time2

    kount_not_found_in_first = 0
    nauxv4 = 0
    call cputim(time1)

    nullify(aux_vec)
    nullify(aux_vec_gat)
    nullify(nauxv4_gat)
    nullify(missing_find)
    nullify(level_found)
    nullify(xcoor)

    fact( 1) = 1.0_rp
    fact( 2) = 1.001_rp
    fact( 3) = 0.999_rp
    fact( 4) = 1.01_rp
    fact( 5) = 0.99_rp
    fact( 6) = 0.5_rp
    fact( 7) = 0.25_rp
    fact( 8) = 0.125_rp
    fact( 9) = 0.0625_rp
    fact(10) = 1.0_rp/32.0_rp
    fact(11) = 1.0_rp/64.0_rp

    call memory_alloca(mem_modul(1:2,modul),'MISSING_FIND','ker_adapel',missing_find,nboun,'DO_NOT_INITIALIZE')
    call memory_alloca(mem_modul(1:2,modul),'LEVEL_FOND'  ,'ker_adapel',level_found ,nboun,'DO_NOT_INITIALIZE')
    do iboun = 1,nboun
       missing_find(iboun) = 1_ip
       level_found(iboun)  = ntest
    end do

    call messages_live('ADAPTIVE WALL EXCHANGE','START SECTION')

    test: do itest = 1,ntest

       !------------------------------------------------------------------
       !
       ! Loop over boundaries - preliminary just to obtain KOUNT
       ! and allocate XCOOR
       !
       !------------------------------------------------------------------

       kount = 0_ip
       do iboun = 1,nboun
          if( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) then
             pblty = ltypb(iboun) 
             kount = kount + ngaus(pblty) 
          end if
       end do

       call memory_alloca(mem_modul(1:2,modul),'XCOOR','ker_adapel',xcoor,ndime,kount)

       !------------------------------------------------------------------
       !
       ! Coordinates XCOOR where to interpolate
       !
       !------------------------------------------------------------------

       kount    = 0_ip
       dist_aux = dexlo_ker
       boundaries: do iboun = 1,nboun
          if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) then
             dist_aux = ywalb(iboun)     ! Variable wall distance (use ywalb instead of dexlo_ker)
             !
             ! Element properties and dimensions
             !
             pblty = ltypb(iboun) 
             pnodb = nnode(pblty)
             ielem = lelbo(iboun)
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             pgaub = ngaus(pblty) 
             pgaus = ngaus(pelty)
             !
             ! Gather operations: ELCOD, BOCOD
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)             
             end do
             do inodb = 1,pnodb 
                ipoin = lnodb(inodb,iboun)
                bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
             end do

             gauss_points: do igaub = 1,pgaub
                !
                ! Obtain normal (baloc(:,ndime) to the surface (following nsi_bouset)
                !                 
                call bouder(&
                     pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                     bocod,baloc,eucta)                                   ! and Jacobian
                call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                kount = kount + 1_ip

                gbcod = 0.0_rp
                do inodb = 1,pnodb
                   gbcod(1:ndime) = gbcod(1:ndime)        &
                        + elmar(pblty)%shape(inodb,igaub) * bocod(1:ndime,inodb)
                end do
                xcoor(1:ndime,kount) = gbcod(1:ndime) - fact(itest) * dist_aux * baloc(1:ndime,ndime) 
             end do gauss_points

          end if

       end do boundaries

       !------------------------------------------------------------------
       !
       ! Interpolation
       !
       !------------------------------------------------------------------

       call messages_live('COMPUTE INTERPOLATION FOR ADAPTIVE STEP '//integer_to_string(itest))

       call fake_waexlo % init       ()
       call fake_waexlo % input      (search_waexlo_seq % method,search_waexlo_par % method,&
            &                        COMM=PAR_COMM_MY_CODE,&
            &                        INTERPOLATION_METHOD=INT_ELEMENT_INTERPOLATION,&
            &                        NAME='ADAPTIVE WALL EXCHANGE '//integer_to_string(itest),&
            &                        DERIVATIVES=.true.)
       call fake_waexlo % preprocess (xcoor,meshe(ndivi))

       !call PAR_BARRIER()
       !print*,'end 1=',kfl_paral,memory_size(fake_waexlo % found),associated(fake_waexlo % matrix)
       !call runend('O.K.!')
       !call fake_waexlo % output     (UNIT=momod(modul) % lun_outpu,HEADER='        ')

       !------------------------------------------------------------------
       !
       ! check if exchange loc has been found for all gauss point on iboun
       !
       !------------------------------------------------------------------

       kount = 0_ip
       boun_check: do iboun = 1,nboun

          if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) then
             !
             ! Element properties and dimensions
             !
             pblty = ltypb(iboun) 
             pnodb = nnode(pblty)
             ielem = lelbo(iboun)
             pgaub = ngaus(pblty) 
             !
             ! Mark the boundary as found - if some of its gauss points is not found it will be changed to 1 in some lines
             !
             missing_find(iboun) = 0  
             do igaub = 1,pgaub
                kount = kount + 1
                if( .not. fake_waexlo % found(kount) ) then
                   if( itest == ntest ) print*,' LAST TEST and point not found with coordinates',xcoor(:,kount)
                   !ACA probablemente mejor escribir la coord de pto de gauss que no encuentra
                   missing_find(iboun) = 1
                end if
             end do
             if ( missing_find(iboun) == 0 ) then  ! If all gauss point in the boundary have been found
                ywalb(iboun)       = ywalb(iboun) * fact(itest)
                level_found(iboun) = itest
                if ( itest > 1 ) kount_not_found_in_first = kount_not_found_in_first + pgaub   ! to be used late for visualization
             end if
          end if

       end do boun_check

       !call PAR_BARRIER()
       !print*,'end 2 =',kfl_paral

       !
       ! Deallocate
       !
       call fake_waexlo % deallo()
       call memory_deallo(mem_modul(1:2,modul),'XCOOR','ker_adapel',xcoor)

       kfl_all_found = 1
       all_found: do iboun = 1,nboun
          if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. missing_find(iboun) == 1 ) kfl_all_found = 0
       end do all_found
       call PAR_MIN(kfl_all_found,'IN MY CODE')
       if( kfl_all_found == 1 ) then
          exit test
       else if( itest == ntest ) then
          call messages_live('ALL POSITIONS HAVE BEEN TESTED AND SOME POINTS HAVE STILL NOT BEEN FOUND')
       end if

    end do test

    !if( IPARALL ) then
    if( 1==1 ) then
       !
       ! Visualization
       !
       kount = 0_ip
       if( INOTMASTER .and. kount_not_found_in_first > 0 ) then
          allocate(aux_vec(n1_aux_vec,kount_not_found_in_first))
       end if
       !call PAR_BARRIER()
       !print*,'end 3 =',kfl_paral

       if( INOTSLAVE ) then
          !
          ! Open file
          !
          fil_tempo = 'exchange_not_found.csv'
          call iofile(0_ip,lun_tempo,fil_tempo,'SYSTEM INFO')
          write(lun_tempo,*)'gp_x,gp_y,gp_z,found_x,found_y,found_z,level,inv_fact,kfl_paral'
          allocate( nauxv4_gat(0:PAR_CODE_SIZE-1) )
       end if
       !call PAR_BARRIER()
       !print*,'end 4 =',kfl_paral
       !
       ! How to use this in Paraview
       ! Open mesh and also open this csv file
       ! Run the filter Filters/ Alphabetical/ Table To Points.   
       ! Tell ParaView what columns are the X, Y and Z coordinate(gp_x,gp_y,gp_z). Be sure to not skip this step. Apply.
       ! Mark keep all data arrays
       ! ParaView probably didn't open up a 3d window (this is a bug). Split screen Horizontal (Icon, top right). 3D View.
       ! You can change the point size so that it looks nicer. 
       ! I believe the best value to plot is inv_fact
       !
       visualization: do iboun = 1,nboun
          if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 .and. level_found(iboun) > 1 ) then
             !
             ! Element properties and dimensions
             !
             pblty = ltypb(iboun) 
             pnodb = nnode(pblty)
             ielem = lelbo(iboun)
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             pgaub = ngaus(pblty) 
             pgaus = ngaus(pelty) 
             !
             ! Gather operations: ELCOD, BOCOD
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin) 
             end do
             do inodb = 1,pnodb 
                ipoin = lnodb(inodb,iboun)
                bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
             end do

             gauss_points_vis: do igaub = 1,pgaub
                !
                ! Obtain normal (baloc(:,ndime) to the surface (following nsi_bouset)
                !
                kount = kount + 1
                call bouder(&
                     pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                     bocod,baloc,eucta)                                   ! and Jacobian
                call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                gbcod = 0.0_rp
                do inodb = 1,pnodb
                   gbcod(1:ndime) = gbcod(1:ndime)        &
                        + elmar(pblty)%shape(inodb,igaub) * bocod(1:ndime,inodb)
                end do

                if( ndime == 3 ) then
                   aux_vec(1:ndime,kount)   = gbcod(1:ndime)
                   aux_vec(4:ndime+3,kount) = gbcod(1:ndime) - ywalb(iboun) * baloc(1:ndime,ndime)
                   aux_vec(7,kount)         = real(level_found(iboun),rp)
                   aux_vec(8,kount)         = 1.0_rp/fact(level_found(iboun))
                   aux_vec(9,kount)         = real(kfl_paral,rp)
                end if
             end do gauss_points_vis

          end if

       end do visualization
       !
       ! Gather all aux_vec
       ! 
       nauxv4 = int(kount_not_found_in_first,4) * int(n1_aux_vec,4)
       call PAR_GATHER(nauxv4,nauxv4_gat,'IN MY CODE')

       !call PAR_BARRIER()
       !print*,'end 5 =',kfl_paral

       if( INOTSLAVE ) then
          nauxv4_tot = 0
          do ipart = 0,PAR_CODE_SIZE-1
             nauxv4_tot = nauxv4_tot + nauxv4_gat(ipart)
          end do
          allocate( aux_vec_gat(nauxv4_tot) )
       end if
       call PAR_GATHERV(aux_vec,aux_vec_gat,nauxv4_gat,'IN MY CODE')

       call cputim(time2) ; time2 = time2 - time1

       call messages_live('MAX FACTOR USED= '//real_to_string(fact(min(ntest,itest)),'(e13.6)'))
       call messages_live('ADAPTIVE WALL EXCHANGE','END SECTION')
       !
       ! Master outputs aux_vec_gat
       !  
       if( INOTSLAVE ) then
          if( ndime == 3 ) then
             do ii = 1,int(nauxv4_tot,ip)/n1_aux_vec
                ! write(lun_tempo,'(9(e14.7,(a)))') ((aux_vec_gat( i + n1_aux_vec*(ii-1) ),','),i=1,9) ! gfortran did not like it this way
                ! but I do not understand why. In any case it accepts it as written below and it works fine
                write(lun_tempo,'(9(e14.7,(a)))') (aux_vec_gat( i + n1_aux_vec*(ii-1) ),',',i=1,9)
             end do
          end if
       end if
       if( associated(nauxv4_gat)  ) deallocate( nauxv4_gat  )
       if( associated(aux_vec_gat) ) deallocate( aux_vec_gat )
       if( associated(aux_vec)     ) deallocate( aux_vec     )

       if( INOTSLAVE ) then
          call iofile_flush_unit(lun_tempo)
          close(lun_tempo)
       end if
    end if
    !
    ! Deallocate
    !
    call memory_deallo(mem_modul(1:2,modul),'MISSING_FIND','ker_adapel',missing_find)
    call memory_deallo(mem_modul(1:2,modul),'LEVEL_FOUND' ,'ker_adapel',level_found)

    !call PAR_BARRIER()
    !print*,'end 6 =',kfl_paral

  end subroutine ker_adapel

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-17
  !> @brief   Get values
  !> @details Get values 
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_waexlo_getval_22(vv_in,vv_out)

    real(rp), intent(in),    pointer :: vv_in (:,:)
    real(rp), intent(inout), pointer :: vv_out(:,:)

    call interp_waexlo % values(vv_in,vv_out)

  end subroutine ker_waexlo_getval_22

  subroutine ker_waexlo_getval_32(vv_in,vv_out)

    real(rp), intent(in),    pointer :: vv_in (:,:,:)
    real(rp), intent(inout), pointer :: vv_out(:,:)
    real(rp),                pointer :: vv_tmp(:,:)

    if( associated(vv_in) ) then
       vv_tmp => vv_in(:,:,1)
    else
       nullify(vv_tmp)
    end if

    call interp_waexlo % values(vv_tmp,vv_out)

  end subroutine ker_waexlo_getval_32

  subroutine ker_waexlo_getval_21(vv_in,vv_out)

    real(rp), intent(in),    pointer :: vv_in (:,:)
    real(rp), intent(inout), pointer :: vv_out(:)
    real(rp),                pointer :: vv_tmp(:)

    if( associated(vv_in) ) then
       vv_tmp => vv_in(:,1)
    else
       nullify(vv_tmp)
    end if

    call interp_waexlo % values(vv_tmp,vv_out)

  end subroutine ker_waexlo_getval_21

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-17
  !> @brief   Get values
  !> @details Get values 
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_waexlo_getder_22(vv_in,vv_out)

    real(rp), intent(in),    pointer :: vv_in (:,:)
    real(rp), intent(inout), pointer :: vv_out(:,:)
    real(rp),                pointer :: vv_tmp(:)

    if( associated(vv_in) ) then
       vv_tmp => vv_in(:,1)
    else
       nullify(vv_tmp)
    end if

    call interp_waexlo % derivatives(vv_tmp,vv_out)

  end subroutine ker_waexlo_getder_22

end module mod_wall_exchange
!> @}
