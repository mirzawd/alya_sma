!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!>----------------------------------------------------------------------
!> @addtogroup Coupling
!> Subroutines for interpolation when using coupling
!> @{
!> @file    mod_interpolation.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Interpolate values at some coordinates
!> @details Interpolate values at some coordinates. 
!>          Manage the color splitting.
!------------------------------------------------------------------------

module mod_interpolation

  use def_kintyp,         only : ip,rp,lg,r1p,r2p,i1p
  use def_spmat,          only : spmat
  use def_master,         only : ISEQUEN,INOTMASTER,kfl_paral
  use def_master,         only : current_code
  use def_domain,         only : ndime,mnode,nbono
  use def_domain,         only : npoin,nelem,nnode
  use def_domain,         only : lnods,ltype,ltopo
  use def_domain,         only : coord,lnnod,lbono
  use def_domain,         only : mnodb,ltopo,ndimb
  use def_domain,         only : elmar,ngaus,lesub
  use def_domain,         only : lnnob
  use def_master,         only : I_AM_IN_SUBD
  use def_master,         only : I_AM_IN_ZONE
  use def_coupli,         only : typ_color_coupling
  use def_coupli,         only : PROJECTION
  use def_coupli,         only : STRESS_PROJECTION
  use def_coupli,         only : BETWEEN_SUBDOMAINS
  use def_coupli,         only : UNKNOWN
  use def_coupli,         only : DIRICHLET_EXPLICIT
  use def_coupli,         only : DIRICHLET_IMPLICIT
  use def_coupli,         only : memor_cou
  use def_coupli,         only : mcoup
  use def_coupli,         only : nboun_cou
  use def_coupli,         only : lnodb_cou 
  use def_coupli,         only : ltypb_cou
  use def_coupli,         only : lnnob_cou
  use def_coupli,         only : lboel_cou
  use def_coupli,         only : lelbo_cou
  use def_coupli,         only : VALUES_ON_NODES
  use def_coupli,         only : VALUES_ON_ELEMENTS
  use def_coupli,         only : VALUES_ON_BOUNDARIES
  use mod_parall,         only : par_bin_part
  use mod_parall,         only : par_bin_size
  use mod_parall,         only : par_bin_boxes
  use mod_parall,         only : par_bin_comin
  use mod_parall,         only : par_bin_comax
  use mod_parall,         only : PAR_COMM_CURRENT 
  use mod_parall,         only : PAR_COMM_COLOR_PERM
  use mod_parall,         only : PAR_MY_CODE_RANK
  use mod_parall,         only : PAR_COMM_COLOR
  use mod_parall,         only : par_part_in_color
  use mod_parall,         only : color_target
  use mod_parall,         only : color_source
  use mod_parall,         only : par_color_to_subd
  use mod_parall,         only : par_color_to_zone
  use mod_parall,         only : par_color_to_code
  use mod_parall,         only : I_AM_IN_COLOR
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_BARRIER
  use mod_communications, only : PAR_START_NON_BLOCKING_COMM
  use mod_communications, only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use mod_communications, only : PAR_END_NON_BLOCKING_COMM
  use mod_integration_rules, only : integration_rules_trapezoidal
  use mod_elmgeo,            only : elmgeo_shapf_deriv_heslo
  use mod_matrix,            only : nullify_spmat
  use mod_matrix,            only : matrix_COO_SpMV_coupling
  use mod_matrix,            only : matrix_COO_SPGEMM
  use def_coupli,            only : ELEMENT_INTERPOLATION
  use def_coupli,            only : BOUNDARY_INTERPOLATION
  use def_coupli,            only : NEAREST_BOUNDARY_NODE
  use def_coupli,            only : NEAREST_ELEMENT_NODE
  use def_coupli,            only : GLOBAL_NUMBERING
  use def_coupli,            only : ON_FLOATING_POINTS
  use def_coupli,            only : FLOATING_TARGET_ENTITY
  use mod_bouder


  implicit none

  character(50), parameter :: vacal = 'mod_interpolation'
  
  private 
  !
  ! Interpolate point values from anywhere
  !
  interface COU_GET_INTERPOLATE_POINTS_VALUES
     module procedure COU_GET_INTERPOLATE_POINTS_VALUES_IP_1,&
          &           COU_GET_INTERPOLATE_POINTS_VALUES_IP_2,&
          &           COU_GET_INTERPOLATE_POINTS_VALUES_RP_0,&
          &           COU_GET_INTERPOLATE_POINTS_VALUES_RP_1,&
          &           COU_GET_INTERPOLATE_POINTS_VALUES_RP_2,&
          &           COU_GET_INTERPOLATE_POINTS_VALUES_RP_3
  end interface COU_GET_INTERPOLATE_POINTS_VALUES

  public :: COU_GET_INTERPOLATE_POINTS_VALUES                 ! Interpolate values afetr initialization
  public :: COU_TARGET_INTERFACE_MASS_MATRIX
  public :: COU_SOURCE_INTERFACE_MASS_MATRIX
  public :: COU_PROJECTION_TYPE
  public :: COU_GENERATE_PROJECTION_MATRIX
  
contains 

  subroutine COU_GET_INTERPOLATE_POINTS_VALUES_IP_1(xx_in,xx_recv,coupling)
    
    integer(ip),              pointer, intent(in)    :: xx_in(:)
    integer(ip),              pointer, intent(inout) :: xx_recv(:)
    type(typ_color_coupling),          intent(inout) :: coupling
    integer(ip)                                      :: ndofn
    integer(ip)                                      :: nsize_in,ii
    integer(ip)                                      :: nsize_recv
    real(rp),                 pointer                :: xxr_in(:)
    real(rp),                 pointer                :: xxr_recv(:)

    nullify(xxr_in)
    nullify(xxr_recv)
    
    ndofn      = 1_ip
    nsize_in   = memory_size(xx_in)
    nsize_recv = memory_size(xx_recv)

    call memory_alloca(memor_cou,'xxr_in'  ,vacal,xxr_in  ,max(1_ip,nsize_in))
    call memory_alloca(memor_cou,'xxr_recv',vacal,xxr_recv,max(1_ip,nsize_recv))

    do ii = 1,nsize_in
       xxr_in(ii) = real(xx_in(ii),rp)
    end do    
    call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xxr_in,xxr_recv,coupling)
     do ii = 1,nsize_recv
       xx_recv(ii) = int(xxr_recv(ii),ip)
    end do

    call memory_deallo(memor_cou,'xxr_in'  ,vacal,xxr_in  )
    call memory_deallo(memor_cou,'xxr_recv',vacal,xxr_recv)
    
  end subroutine COU_GET_INTERPOLATE_POINTS_VALUES_IP_1

  subroutine COU_GET_INTERPOLATE_POINTS_VALUES_IP_2(xx_in,xx_recv,coupling)
    
    integer(ip),              pointer, intent(in)    :: xx_in(:,:)
    integer(ip),              pointer, intent(inout) :: xx_recv(:,:)
    type(typ_color_coupling),          intent(inout) :: coupling
    integer(ip)                                      :: ndofn
    integer(ip)                                      :: nsize_in,ii
    integer(ip)                                      :: nsize_recv
    real(rp),                 pointer                :: xxr_in(:,:)
    real(rp),                 pointer                :: xxr_recv(:,:)

    nullify(xxr_in)
    nullify(xxr_recv)
    
    ndofn      = max(memory_size(xx_in,1_ip),memory_size(xx_recv,1_ip))
    nsize_in   = memory_size(xx_in,2_ip)
    nsize_recv = memory_size(xx_recv,2_ip)

    call memory_alloca(memor_cou,'xxr_in'  ,vacal,xxr_in  ,max(1_ip,ndofn),max(1_ip,nsize_in))
    call memory_alloca(memor_cou,'xxr_recv',vacal,xxr_recv,max(1_ip,ndofn),max(1_ip,nsize_recv))

    do ii = 1,nsize_in
       xxr_in(:,ii) = real(xx_in(:,ii),rp)
    end do    
    call COU_GET_INTERPOLATE_POINTS_VALUES_RP(max(1_ip,ndofn),xxr_in,xxr_recv,coupling)
     do ii = 1,nsize_recv
       xx_recv(:,ii) = int(xxr_recv(:,ii),ip)
    end do

    call memory_deallo(memor_cou,'xxr_in'  ,vacal,xxr_in  )
    call memory_deallo(memor_cou,'xxr_recv',vacal,xxr_recv)
    
  end subroutine COU_GET_INTERPOLATE_POINTS_VALUES_IP_2

  subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_0(ndofn,xx_in,xx_recv,coupling)
    implicit none
    integer(ip),                       intent(in)    :: ndofn
    real(rp),                          intent(in)    :: xx_in(ndofn,*)
    real(rp),                          intent(inout) :: xx_recv(ndofn,*)
    type(typ_color_coupling),          intent(inout) :: coupling

    call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,xx_recv,coupling)       

  end subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_0

  subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_1(xx_in,xx_recv,coupling)
    implicit none
    real(rp),                 pointer, intent(in)    :: xx_in(:)
    real(rp),                 pointer, intent(inout) :: xx_recv(:)
    type(typ_color_coupling),          intent(inout) :: coupling
    real(rp)                                         :: yy_in(2)
    real(rp)                                         :: yy_recv(2)
    integer(ip)                                      :: ndofn
    ndofn = 1_ip

    !if( associated(xx_in) ) then
    !   if(  size(xx_in) < npoin .and. &
    !        coupling % kfl_source_value == VALUES_ON_NODES .and. &
    !        PAR_MY_CODE_RANK /= 0 ) call runend('WRONG SIZE 1: CANNOT INTERPOLATE')
    !end if

    if( .not. associated(xx_in) .and. associated(xx_recv) ) then
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,yy_in,xx_recv,coupling)
    else if( .not. associated(xx_in) .and. .not. associated(xx_recv) ) then
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,yy_in,yy_recv,coupling)
    else if( associated(xx_in) .and. .not. associated(xx_recv) ) then
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,yy_recv,coupling)
    else 
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,xx_recv,coupling)
    end if
  end subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_1

  subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_2(xx_in,xx_recv,coupling)

    real(rp),                 pointer, intent(in)    :: xx_in(:,:)
    real(rp),                 pointer, intent(inout) :: xx_recv(:,:)
    type(typ_color_coupling),          intent(inout) :: coupling
    real(rp)                                         :: yy_in(2)
    real(rp)                                         :: yy_recv(2)
    integer(ip)                                      :: ndofn

    !if( associated(xx_in) ) then
    !   if(  size(xx_in,2) < npoin .and. &
    !        coupling % kfl_source_value == VALUES_ON_NODES .and. &
    !        PAR_MY_CODE_RANK /= 0 ) call runend('WRONG SIZE 2: CANNOT INTERPOLATE')
    !end if

    if( .not. associated(xx_in) .and. associated(xx_recv) ) then
       ndofn = size(xx_recv,1)
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,yy_in,xx_recv,coupling)
    else if( .not. associated(xx_in) .and. .not. associated(xx_recv) ) then
       ndofn = 1
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,yy_in,yy_recv,coupling)
    else if( associated(xx_in) .and. .not. associated(xx_recv) ) then
       ndofn = size(xx_in,1)
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,yy_recv,coupling)
    else 
       ndofn = max(memory_size(xx_in,1_ip),memory_size(xx_recv,1_ip))
       if( coupling % commd % lrecv_dim > 0 .and. coupling % commd % lsend_dim > 0 ) then
          if( ndofn /= size(xx_recv,1) .and. PAR_MY_CODE_RANK /= 0 ) call runend('WRONG SIZE 3: CANNOT INTERPOLATE')
       end if
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,xx_recv,coupling)
    end if

  end subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_2

  subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_3(xx_in,xx_recv,coupling)
    implicit none
    real(rp),                 pointer, intent(in)    :: xx_in(:,:,:)
    real(rp),                 pointer, intent(inout) :: xx_recv(:,:)
    type(typ_color_coupling),          intent(inout) :: coupling
    real(rp)                                         :: yy_in(2)
    real(rp)                                         :: yy_recv(2)
    integer(ip)                                      :: ndofn

    !if( associated(xx_in) ) then
    !   if( size(xx_in,2) < npoin .and. PAR_MY_CODE_RANK /= 0 ) then
    !      call runend('WRONG SIZE 4: CANNOT INTERPOLATE')
    !   end if
    !end if

    if( .not. associated(xx_in) .and. associated(xx_recv) ) then
       ndofn = size(xx_recv,1)
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,yy_in,xx_recv,coupling)
    else if( .not. associated(xx_in) .and. .not. associated(xx_recv) ) then
       ndofn = 1
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,yy_in,yy_recv,coupling)
    else if( associated(xx_in) .and. .not. associated(xx_recv) ) then
       ndofn = size(xx_in,1)
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,yy_recv,coupling)
    else 
       ndofn = max(memory_size(xx_in,1_ip),memory_size(xx_recv,1_ip))
       if( coupling % commd % lrecv_dim > 0 .and. coupling % commd % lsend_dim > 0 ) then
          if( ndofn /= size(xx_recv,1) .and. PAR_MY_CODE_RANK /= 0 ) call runend('WRONG SIZE 5: CANNOT INTERPOLATE')
       end if
       call COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,xx_recv,coupling)
    end if

  end subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP_3


  subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP(ndofn,xx_in,xx_out,coupling,what_to_do)

    implicit none
    integer(ip),              intent(in)           :: ndofn
    real(rp),                 intent(in)           :: xx_in(ndofn,*)
    real(rp),                 intent(out)          :: xx_out(ndofn,*)
    type(typ_color_coupling), intent(inout)        :: coupling
    character(*),             intent(in), optional :: what_to_do
    integer(ip)                                    :: ineig
    real(rp),    pointer                           :: xx_recv(:,:)
    real(rp)                                       :: time1,time2,time3,weight
    logical(lg)                                    :: done_ini
    integer(ip)                                    :: ipoin, idofn,kpoin
    integer(ip)                                    :: color_target_aux
    character(100), PARAMETER                      :: vacal = "COU_GET_INTERPOLATE_POINTS_VALUES"
    !
    ! Allocate arrays
    !
    nullify(xx_recv)
    call memory_alloca(memor_cou,'xx_recv',vacal,xx_recv,ndofn,max(1_ip,coupling % commd % lrecv_dim))
    color_target_aux = color_target
    !
    ! Send and receive 
    !
    call cputim(time1)
    call PAR_SEND_RECEIVE_TO_ALL(ndofn,xx_in,xx_recv,coupling % commd,'ASYNCHRONOUS',PERMUTE_SEND=.true.)!,PERMUTE_RECV=.true.)
    call cputim(time2)
    !  
    ! Operate Transmission matrix
    !
    if(associated(coupling % ltransmat_target)) then
       done_ini = .false.
       do ineig = 1,coupling % commd % nneig
          if(associated(coupling % ltransmat_target(ineig) % iA)) then
                if(.not. done_ini)   xx_out(1:ndofn,1:coupling % ltransmat_target(ineig) % nrows ) = 0.0_rp
                call matrix_COO_SpMV_coupling(coupling % ltransmat_target(ineig),ndofn,xx_recv,xx_out,noini=.true.)
                done_ini = .true.
          end if
       end do
    end if
    call cputim(time3)
    !
    ! Deallocate arrays
    !
    call memory_deallo(memor_cou,'xx_recv',vacal,xx_recv)
    
    if(  coupling % kfl_par_transm == 0                      .and.  &
       & current_code == coupling % code_target              .and.  &
       & coupling  % target_entity  /= FLOATING_TARGET_ENTITY )  then


       call memory_alloca(memor_cou,'xx_recv',vacal,xx_recv,ndofn,max(1_ip,npoin))
     
       do kpoin = 1, coupling % wet % npoin_wet

          ipoin  = coupling % wet % lpoin_wet(kpoin)

          if(  coupling % itype == ELEMENT_INTERPOLATION  .or. &
               coupling % itype == NEAREST_BOUNDARY_NODE  .or. &
               coupling % itype == BOUNDARY_INTERPOLATION .or. &
               coupling % itype == NEAREST_ELEMENT_NODE   .or. &
               coupling % itype == GLOBAL_NUMBERING       ) then

             weight = coupling % wet % weight_wet_imp(kpoin)
          else
             weight = 1.0_rp
          endif
          do idofn = 1,ndofn
             xx_recv(idofn,ipoin) = xx_out(idofn,kpoin)*weight
          end do
       end do
       color_target = coupling  % color_target
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,xx_recv,'SUM','IN CURRENT TARGET COLOR')
       do kpoin = 1, coupling % wet % npoin_wet
          ipoin  = coupling % wet % lpoin_wet(kpoin)
          do idofn = 1,ndofn
             xx_out(idofn,kpoin) = xx_recv(idofn,ipoin)
          end do
       end do
       call memory_deallo(memor_cou,'xx_recv',vacal,xx_recv)
    endif
    ! 
    ! Timers
    !
    coupling % cputim(8) = coupling % cputim(8) + time2 - time1
    coupling % cputim(9) = coupling % cputim(9) + time3 - time2
    color_target = color_target_aux

  end subroutine COU_GET_INTERPOLATE_POINTS_VALUES_RP
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Residual projection
  !> @details Residual projection. There are two possibilities.
  !>          F = nodal force
  !>          f = stress = F/S
  !>          
  !>          1. Find nodal stress on source fi,s = Fi,s / Mii
  !>
  !>          2. Compute force Fi = 
  !>                            +-
  !>                      Fi =  |  ( Ni,s f_source ) ds
  !>                           -+ source/target
  !>                            +-
  !>             Option 1:   =  |  ( Ni,t f_source ) ds  (f_source interpolated on target)
  !>                           -+ target
  !>                            +-
  !>             Option 2:   =  |  ( Ni,t Nj,s fj,s ) ds (Ni interpolated on source)
  !>                           -+ source
  !>
  !>          \verbatim
  !>
  !>          1. All on target
  !>          ----------------
  !>
  !>             1. On source: compute nodal stress fj_source = Fj_s / Mjj 
  !>             2. On source: interpolate f_source at target Gauss points
  !>             3. On target:
  !>             +-                   +-
  !>             | ( Ni Nj fj ) ds =  |  ( Ni f_source ) ds
  !>            -+ target            -+ target 
  !>                                  +-                      +--  +--  +--
  !>                            Fi =  |  ( Ni f_source ) ds = \    \    \   wg |J| Ni f_source
  !>                                 -+ source                /__  /__  /__
  !>                                                         iboun  i    g
  !>
  !>            NB: If a closed rule is selected and mesh coincide
  !>            
  !>                  +-                      
  !>            Fi =  |  ( Ni Nj Fj,s / Mjj ) ds = Ni Ni Fi,s / Mii * Mii = Fi,s
  !>                 -+ source            
  !>          
  !>
  !>             o-------o-------o-------o-------o source
  !>                  ||    ||        ||    ||  
  !>                  \/    \/        \/    \/  
  !>             x----g-----g----x----g-----g----x target
  !>
  !>              <--- iboun ---> <--- iboun --->
  !>                          
  !>          2. On source and target
  !>          -----------------------
  !>
  !>             +-                   +-
  !>             | ( Ni Nj fj ) ds =  |  ( Ni f_source ) ds
  !>            -+ target            -+ source 
  !>                                  +-
  !>                        Mii fi =  |  ( Ni F_source / M_source ) ds
  !>                                 -+ source 
  !>
  !>
  !>             o-g---g-o-------o-------o-------o source <= RHS computed here
  !>               /\
  !>             i || Ni
  !>             x---------------x---------------x target <= Mii computed here
  !>
  !>             The problem with this option is that in the latter
  !>             example all the CPUs of the source would require
  !>             the boundaries of the target as they all have
  !>             influence on the middle node
  !>
  !>             In this case, wet points are Gauss points
  !>
  !>          -------------------
  !>          1. Compute residual
  !>          -------------------
  !>
  !>                    CPU3                  CPU4
  !>          RA                  RB3' RB4'               RC
  !>          o---------------------o o-------------------o  source
  !>
  !>
  !>          --------------------
  !>          2. Exchange residual: Rb=Rb3+Rb4 
  !>          --------------------
  !>
  !>                    CPU3                  CPU4
  !>          RA                  RB3 RB4                 RC
  !>          o---------------------o o-------------------o  source
  !>          <=====================>
  !>                   SAB
  !>
  !>
  !>          ---------------------
  !>          3. Interpolate stress (at target Gauss points x)
  !>          --------------------- 
  !>
  !>                    CPU3                  CPU4
  !>
  !>          o----x------x---------o ox-----x-----------o  source
  !>           I(R)/SAB I(R)/SAB   I(R)/SBC I(R)/SBC                          
  !>
  !>
  !>          ------------------
  !>          4. Assemble stress (integrate on target)
  !>          ------------------ 
  !>                       
  !>                 CPU1               CPU2
  !>          ra            rb1' rb2'            rc
  !>          o----x------x----o o----x------x----o         target
  !>
  !>
  !>          ------------------
  !>          5. Exchange stress: rb = rb1' + rb2'
  !>          ------------------ 
  !>                       
  !>                 CPU1               CPU2
  !>          ra              rb rb            rc
  !>          o----x------x----o o----x------x----o         target
  !>
  !>          \endverbatim
  !>
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Assemble Gauss point contributions
  !> @details 
  !>           X_RECV_TMP
  !>              ||
  !>              \/        XX_OUT(NPOIN_WET)
  !>          o---g-----g---o
  !>
  !>
  !----------------------------------------------------------------------

  subroutine COU_PROJECTION(itask,coupling,ndofn,xx_recv_tmp,xx_out)
    integer(ip),              intent(in)  :: itask
    type(typ_color_coupling), intent(in)  :: coupling
    integer(ip),              intent(in)  :: ndofn
    real(rp),        pointer, intent(in)  :: xx_recv_tmp(:,:)
    real(rp),                 intent(out) :: xx_out(ndofn,*)
    integer(ip)                           :: iboun_global,iboun_wet
    integer(ip)                           :: pblty,pnodb,pgaub,igaub
    integer(ip)                           :: inodb,kgaub,jgaub
    integer(ip)                           :: ipoin_wet,npoin_wet
    real(rp)                              :: gbsur

!!$    real(rp), pointer :: xx(:,:)
!!$    if( coupling % geome % npoin_wet > 0 ) then
!!$       !
!!$       ! Initalize residual
!!$       !
!!$       allocate(xx(ndofn,npoin))
!!$       npoin_wet = coupling % geome % npoin_wet
!!$       do ipoin = 1,npoin
!!$          xx(1:ndofn,ipoin) = 0.0_rp
!!$       end do
!!$       jgaub = 0      
!!$       do iboun_wet = 1,coupling % geome % nboun_wet          
!!$          iboun_global = coupling % geome % lboun_wet(iboun_wet)
!!$          pblty        = abs(ltypb_cou(iboun_global))
!!$          pnodb        = lnnob_cou(iboun_global)
!!$          pgaub        = coupling % wet % proje_target(iboun_wet) % pgaub
!!$          do igaub = 1,pgaub
!!$             jgaub = jgaub + 1
!!$             kgaub = coupling % geome % status(jgaub) ! Permutation of wet point             
!!$             gbsur = coupling % wet % proje_target(iboun_wet) % gbsur(igaub)
!!$             do inodb = 1,pnodb
!!$                ipoin = lnodb_cou(inodb,iboun_global)
!!$                ipoin_wet = coupling % wet % proje_target(iboun_wet) % permu(inodb)
!!$                xx(1:ndofn,ipoin) = xx(1:ndofn,ipoin) &
!!$                     & + gbsur * xx_recv_tmp(1:ndofn,kgaub) &
!!$                     & * coupling % wet % proje_target(iboun_wet) % shapb(inodb,igaub)
!!$             end do
!!$          end do          
!!$       end do       
!!$       if( current_code == 2 ) then
!!$          call PAR_INTERFACE_NODE_EXCHANGE(xx,'SUM','IN MY CODE')
!!$       end if       
!!$       do ipoin_wet = 1,npoin_wet
!!$          ipoin = coupling % geome % lpoin_wet(ipoin_wet)
!!$          xx_out(1:ndofn,ipoin_wet) = xx(1:ndofn,ipoin)
!!$       end do
!!$       deallocate( xx )
!!$       if( itask == 2 ) then
!!$          do ipoin_wet = 1,npoin_wet
!!$             xx_out(1:ndofn,ipoin_wet) = xx_out(1:ndofn,ipoin_wet) / coupling % geome % vmasb_wet(ipoin_wet)
!!$          end do
!!$       end if
!!$    end if

    if( coupling % wet % npoin_wet > 0 ) then
       !
       ! Initalize residual
       !
       npoin_wet = coupling % wet % npoin_wet
       do ipoin_wet = 1,npoin_wet
          xx_out(1:ndofn,ipoin_wet) = 0.0_rp
       end do

       jgaub = 0

       do iboun_wet = 1,coupling % wet % nboun_wet

          iboun_global = coupling % wet % lboun_wet(iboun_wet)
          pblty        = abs(ltypb_cou(iboun_global))
          pnodb        = lnnob_cou(iboun_global)
          pgaub        = coupling % wet % proje_target(iboun_wet) % pgaub
          do igaub = 1,pgaub
             jgaub = jgaub + 1
             kgaub = coupling % geome % status(jgaub) ! Permutation of wet point             
             gbsur = coupling % wet % proje_target(iboun_wet) % gbsur(igaub)
             do inodb = 1,pnodb
                ipoin_wet = coupling % wet % proje_target(iboun_wet) % permu(inodb)
                xx_out(1:ndofn,ipoin_wet) = xx_out(1:ndofn,ipoin_wet) &
                     & + gbsur * xx_recv_tmp(1:ndofn,kgaub) &
                     & * coupling % wet % proje_target(iboun_wet) % shapb(inodb,igaub)
             end do
          end do

       end do

       if( itask == 2 ) then
          do ipoin_wet = 1,npoin_wet
             xx_out(1:ndofn,ipoin_wet) = xx_out(1:ndofn,ipoin_wet) / coupling % wet % vmasb_wet(ipoin_wet)
          end do
       end if

    end if

  end subroutine COU_PROJECTION

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Compute mass matrix for this coupling
  !> @details Compute mass matrix for this coupling
  !>
  !----------------------------------------------------------------------

  subroutine COU_TARGET_INTERFACE_MASS_MATRIX(coupling)
    type(typ_color_coupling), intent(inout) :: coupling
    integer(ip)                             :: iboun,kboun,nboun_wet
    integer(ip)                             :: pblty,pnodb,pgaub,igaub
    integer(ip)                             :: inodb,ipoin,npoin_wet
    integer(ip)                             :: ipoin_wet
    real(rp)                                :: baloc(3,3),gbsur,gbdet
    real(rp)                                :: bocod(ndime,mnodb)
    real(rp),                 pointer       :: vmasb(:)

    nboun_wet = coupling % wet % nboun_wet
    npoin_wet = coupling % wet % npoin_wet

    if( associated(coupling % wet % vmasb_wet) ) then
       call memory_deallo(memor_cou,'VMASB','cou_interface_mass_matrix',&
            coupling % wet % vmasb_wet)
    end if
    call memory_alloca(memor_cou,'VMASB','cou_interface_mass_matrix',&
         coupling % wet % vmasb_wet,coupling % wet % npoin_wet)

    allocate( vmasb(npoin) )
    do ipoin = 1,npoin
       vmasb(ipoin) = 0.0_rp
    end do

    do kboun = 1,nboun_wet

       iboun = coupling % wet % lboun_wet(kboun)
       pblty = ltypb_cou(iboun)

       if( pblty > 0 ) then

          pnodb = lnnob_cou(iboun)
          pgaub = ngaus(pblty)
          do inodb = 1,pnodb
             ipoin = lnodb_cou(inodb,iboun)
             bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
          end do

          do igaub = 1,pgaub
             call bouder(pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),bocod,baloc,gbdet)
             gbsur = elmar(pblty) % weigp(igaub) * gbdet 
             do inodb = 1,pnodb
                ipoin = lnodb_cou(inodb,iboun)
                vmasb(ipoin) = vmasb(ipoin) + gbsur * elmar(pblty) % shape(inodb,igaub)
                !kpoin = list_wet_nodes(ipoin)
                !coupling % geome % vmasb_wet(kpoin) = coupling % geome % vmasb_wet(kpoin) &
                !     + gbsur * elmar(pblty) % shape(inodb,igaub)
             end do
          end do

       end if

    end do
    !
    ! OJO: mejorar para que el intercambio se haga solamente en este color
    ! y que sea del tamano npoin-wet...mmmm
    !
    call PAR_INTERFACE_NODE_EXCHANGE(vmasb,'SUM','IN MY CODE')
    do ipoin_wet = 1,coupling % wet % npoin_wet
       ipoin = coupling % wet % lpoin_wet(ipoin_wet)
       coupling % wet % vmasb_wet(ipoin_wet) = vmasb(ipoin) 
    end do

    deallocate(vmasb)

  end subroutine COU_TARGET_INTERFACE_MASS_MATRIX

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Compute mass matrix for this coupling
  !> @details Compute mass matrix for this coupling
  !>
  !----------------------------------------------------------------------

  subroutine COU_SOURCE_INTERFACE_MASS_MATRIX(coupling,mirror_coupling)
    type(typ_color_coupling), intent(inout) :: coupling
    type(typ_color_coupling), intent(inout) :: mirror_coupling
    integer(ip)                             :: iboun
    integer(ip)                             :: pblty,pnodb,pgaub,igaub
    integer(ip)                             :: inodb,ipoin
    integer(ip)                             :: iboun_wet
    real(rp)                                :: baloc(3,3),gbsur,gbdet
    real(rp)                                :: bocod(ndime,mnodb)
    !
    ! Allocate memory
    !
    if( associated(coupling % geome % vmasb) ) then
       call memory_deallo(memor_cou,'VMASB','cou_interface_mass_matrix',coupling % geome % vmasb)
    end if
    call memory_alloca(memor_cou,'VMASB','cou_interface_mass_matrix',coupling % geome % vmasb,npoin) 

    do iboun_wet = 1,mirror_coupling % wet % nboun_wet
       iboun = mirror_coupling % wet % lboun_wet(iboun_wet) 
       pblty = abs(ltypb_cou(iboun))
       pgaub = ngaus(pblty) 
       pnodb = lnnob_cou(iboun)
       do inodb = 1,pnodb
          ipoin = lnodb_cou(inodb,iboun)
          bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
       end do
       do igaub = 1,pgaub
          call bouder(pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),bocod,baloc,gbdet)                   
          gbsur = elmar(pblty) % weigp(igaub) * gbdet 
          do inodb = 1,pnodb
             ipoin = lnodb_cou(inodb,iboun)
             coupling % geome % vmasb(ipoin) = coupling % geome % vmasb(ipoin) &
                  + gbsur * elmar(pblty) % shape(inodb,igaub) 
          end do
       end do
    end do
    !
    ! OJO: mejorar para que el intercambio se haga solamente en este color
    ! y que sea del tamano npoin-wet...mmmm
    !
    call PAR_INTERFACE_NODE_EXCHANGE(coupling % geome % vmasb,'SUM','IN MY CODE')

  end subroutine COU_SOURCE_INTERFACE_MASS_MATRIX

  subroutine COU_PROJECTION_TYPE(coupling)
    type(typ_color_coupling), intent(inout)    :: coupling
    integer(ip)                                :: nboun_wet
    integer(ip)                                :: iboun_wet,iboun_global,pblty
    integer(ip), allocatable                   :: inv_perm(:) 
    integer(ip)                                :: ipoin_wet,ipoin_global,pgaub
    integer(ip)                                :: inodb,pnodb,igaub
    real(rp)                                   :: gbsur,baloc(3,3),bocod(ndime,mnodb)
    real(rp)                                   :: shapb(mnodb),derib(ndimb,mnodb)
    real(rp)                                   :: gbdet
    real(rp),    allocatable                   :: posgp(:,:),weigp(:)
    !
    ! INV_PERM: temporary permutation array
    !
    allocate( inv_perm(npoin) )
    do ipoin_global = 1,npoin
       inv_perm(ipoin_global) = 0
    end do
    do ipoin_wet = 1,coupling % wet % npoin_wet
       ipoin_global = coupling % wet % lpoin_wet(ipoin_wet)
       inv_perm(ipoin_global) = ipoin_wet
    end do
    !
    ! Compute projection type for each wet boundary
    !
    nboun_wet = coupling % wet % nboun_wet
    allocate(coupling % wet % proje_target(nboun_wet) )

    if( coupling % ngaus == 0 ) then
       !
       ! Use kernel default number of Gauss points
       !
       do iboun_wet = 1,nboun_wet
          iboun_global = coupling % wet % lboun_wet(iboun_wet)
          pblty        = abs(ltypb_cou(iboun_global))
          pgaub        = ngaus(pblty)
          pnodb        = nnode(pblty)

          coupling % wet % proje_target(iboun_wet) % pgaub = pgaub
          allocate(coupling % wet % proje_target(iboun_wet) % shapb(pnodb,pgaub) )
          allocate(coupling % wet % proje_target(iboun_wet) % gbsur(pgaub) )
          allocate(coupling % wet % proje_target(iboun_wet) % permu(pnodb) )

          do inodb = 1,pnodb
             ipoin_global = lnodb_cou(inodb,iboun_global)
             bocod(1:ndime,inodb) = coord(1:ndime,ipoin_global)
          end do
          do igaub = 1,pgaub
             call bouder(pnodb,ndime,ndimb,elmar(pblty) % deriv(:,:,igaub),bocod,baloc,gbdet)                
             gbsur = elmar(pblty) % weigp(igaub) * gbdet 
             coupling % wet % proje_target(iboun_wet) % gbsur(igaub) = gbsur
             do inodb = 1,pnodb
                ipoin_global = lnodb_cou(inodb,iboun_global)
                if( ipoin_global == 0 ) call runend('COU_DEFINE_WET_GEOMETRY: WE ARE IN TROUBLE')
                ipoin_wet = inv_perm(ipoin_global)
                coupling % wet % proje_target(iboun_wet) % permu(inodb) = ipoin_wet
                coupling % wet % proje_target(iboun_wet) % shapb(inodb,igaub) = elmar(pblty) % shape(inodb,igaub)
             end do
          end do
       end do
    else
       !
       ! Use a user-prescribed number of Gauss points
       !
       pgaub = coupling % ngaus
       allocate(posgp(ndimb,pgaub),weigp(pgaub)) 
       call integration_rules_trapezoidal(ndimb,pgaub,posgp,weigp)
       do iboun_wet = 1,nboun_wet
          iboun_global = coupling % wet % lboun_wet(iboun_wet)
          pblty        = abs(ltypb_cou(iboun_global))
          pnodb        = nnode(pblty)
          do inodb = 1,pnodb
             ipoin_global = lnodb_cou(inodb,iboun_global)
             bocod(1:ndime,inodb) = coord(1:ndime,ipoin_global)
          end do

          coupling % wet % proje_target(iboun_wet) % pgaub = pgaub
          allocate(coupling % wet % proje_target(iboun_wet) % shapb(pnodb,pgaub) )
          allocate(coupling % wet % proje_target(iboun_wet) % gbsur(pgaub) )
          allocate(coupling % wet % proje_target(iboun_wet) % permu(pnodb) )

          do igaub = 1,pgaub
             call elmgeo_shapf_deriv_heslo(ndimb,pnodb,[posgp(1,igaub)],shapb,derib)
             call bouder(pnodb,ndime,ndimb,derib,bocod,baloc,gbdet)         
             gbsur = weigp(igaub) * gbdet 
             coupling % wet % proje_target(iboun_wet) % gbsur(igaub) = gbsur
             do inodb = 1,pnodb
                ipoin_global = lnodb_cou(inodb,iboun_global)
                if( ipoin_global == 0 ) call runend('COU_DEFINE_WET_GEOMETRY: WE ARE IN TROUBLE')
                ipoin_wet = inv_perm(ipoin_global)
                coupling % wet % proje_target(iboun_wet) % permu(inodb) = ipoin_wet
                coupling % wet % proje_target(iboun_wet) % shapb(inodb,igaub) = shapb(inodb)
             end do
          end do

       end do
       deallocate(posgp,weigp)
    end if

    deallocate(inv_perm)

  end subroutine COU_PROJECTION_TYPE

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/07/2016
  !> @brief   Generate matrix to assemble Gauss point contributions
  !>
  !----------------------------------------------------------------------

  subroutine COU_GENERATE_PROJECTION_MATRIX(coupling,M)

    implicit none
    type(typ_color_coupling), intent(in)  :: coupling
    type(spmat), intent(inout)            :: M
    integer(ip)                           :: itask
    integer(ip)                           :: iboun_global,iboun_wet
    integer(ip)                           :: pnodb,pgaub,igaub
    integer(ip)                           :: inodb,kgaub,jgaub
    integer(ip)                           :: ipoin_wet
    real(rp)                              :: gbsur, val
    integer(ip)                           :: nentr

    real(rp),PARAMETER :: tolval = 1.0e-10_rp

    if( coupling % itype == STRESS_PROJECTION ) then
       itask = 1
    else if( coupling % itype == PROJECTION ) then
       itask = 2
    else
       call runend("COUPLING MUST BE A PROJECTION")
    endif

    if( coupling % wet % npoin_wet > 0 ) then
       !
       ! Count entries in M
       !
       nentr = 0
       jgaub = 0
       do iboun_wet = 1,coupling % wet % nboun_wet
          iboun_global = coupling % wet % lboun_wet(iboun_wet)
          pnodb        = lnnob_cou(iboun_global)
          pgaub        = coupling % wet % proje_target(iboun_wet) % pgaub
          do igaub = 1,pgaub
             jgaub = jgaub + 1
             gbsur = coupling % wet % proje_target(iboun_wet) % gbsur(igaub)
             do inodb = 1,pnodb
                ipoin_wet = coupling % wet % proje_target(iboun_wet) % permu(inodb)
                val = gbsur * coupling % wet % proje_target(iboun_wet) % shapb(inodb,igaub) 
                if(abs(val) > tolval) then
                   nentr = nentr + 1
                endif
             enddo
          enddo
       end do
       !
       ! Allocate M
       !
       call nullify_spmat(M)
       call memory_alloca(memor_cou,'M','COU_GENERATE_PROJECTION_MATRIX',M, 1_ip, nentr)
       !
       ! Evaluate coeficients of M
       !
       M % ndof1 = 1
       M % ndof2 = 1
       M % nrows = coupling % wet % npoin_wet
       M % ncols = coupling % wet % number_wet_points 
       jgaub = 0
       nentr = 1
       do iboun_wet = 1,coupling % wet % nboun_wet

          iboun_global = coupling % wet % lboun_wet(iboun_wet)
          pnodb        = lnnob_cou(iboun_global)
          pgaub        = coupling % wet % proje_target(iboun_wet) % pgaub
          do igaub = 1,pgaub
             jgaub = jgaub + 1
             kgaub = coupling % geome % status(jgaub) ! Permutation of wet point             
             gbsur = coupling % wet % proje_target(iboun_wet) % gbsur(igaub)
             do inodb = 1,pnodb
                ipoin_wet = coupling % wet % proje_target(iboun_wet) % permu(inodb)
                val = gbsur * coupling % wet % proje_target(iboun_wet) % shapb(inodb,igaub)
                if(abs(val) > tolval) then
                   if(itask == 2) val = val/coupling % wet % vmasb_wet(ipoin_wet)
                   M % iA(nentr)      = ipoin_wet
                   M % jA(nentr)      = kgaub
                   M % vA(1,1,nentr)  = val 
                   nentr = nentr + 1
                endif
             end do
          end do
       end do
    end if
  
 end subroutine COU_GENERATE_PROJECTION_MATRIX

end module mod_interpolation
!> @}
