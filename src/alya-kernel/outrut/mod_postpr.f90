!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @name    Toolbox for pstprocess
!> @file    mod_postpr.f90
!> @author  houzeaux
!> @date    2018-02-19
!> @brief   postprocess
!> @details postprocess
!>          There are two formats:
!>          Alya binary: *.alyabin
!>          MPIO binary format: *.mpio.bin
!> @{
!-----------------------------------------------------------------------

module mod_postpr

  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use def_postpr
  use mod_mpio_seq_log
  use mod_mpio_postpr
  use mod_permut
  use def_inpout
  use def_mpio
  use mod_postpr_tools
  use def_kintyp_mesh,           only : mesh_type_basic
  use mod_outfor,                only : outfor
  use mod_messages,              only : livinf
  use mod_iofile,                only : iofile_flush_unit
  use mod_iofile,                only : iofile_file_exists
  use mod_memory,                only : memory_alloca
  use mod_memory,                only : memory_deallo
  use mod_memory,                only : memory_size
  use mod_communications,        only : PAR_SEND_RECEIVE
  use mod_communications,        only : PAR_SUM
  use mod_communications,        only : PAR_MAX
  use mod_communications_global, only : PAR_BROADCAST
  use mod_opfpos,                only : opfpos_name
  use mod_messages,              only : messages_live
  use mod_postpr_on_basic_mesh,  only : postpr_mesh_mpio
  use mod_postpr_on_basic_mesh,  only : postpr_mesh_alya
  use mod_interp_fe,             only : interp_fe
  use mod_interp_fe,             only : interp_fe_deallocate
  implicit none

  private

  integer(ip)             :: nenti

  character(5)            :: bench_rw
  integer(ip)             :: kfl_perme
  integer(ip)             :: kfl_permb

  integer(ip)             :: bridge_ip(1)
  integer(ip), pointer    :: bridge_ip_vector(:,:)
  integer(ip), pointer    :: bridge_ip_scalar(:)
  real(rp)                :: bridge_rp(1)
  real(rp),    pointer    :: bridge_rp_vector(:,:)
  real(rp),    pointer    :: bridge_rp_scalar(:)

  integer(ip)             :: id_loc
  integer(ip)             :: modul_loc
  integer(ip)             :: enti_posit_loc
  integer(ip)             :: comp_posit_loc
  integer(ip)             :: time_posit_loc
  character(5)            :: wopos_loc(5)
  !
  ! Interpolation variables
  !
  real(rp),    pointer    :: xx_sca(:)
  real(rp),    pointer    :: xx_vec(:,:)
  
  interface postpr
     module procedure &
          postpr_real_scalar , postpr_real_vector , &
          postpr_int_scalar  , postpr_int_vector  , &
          posr3p,posvex
  end interface postpr

  interface postpr_postprocess
     module procedure &
          postpr_postprocess_RP_1 , postpr_postprocess_RP_2 , postpr_postprocess_RP_3 , &
          postpr_postprocess_IP_1 , postpr_postprocess_IP_2 , postpr_postprocess_IP_3 , &
          postpr_postprocess_R3P_1
  end interface postpr_postprocess

  interface postpr_int_to_real
     module procedure &
          postpr_int_to_real_1,postpr_int_to_real_2,postpr_int_to_real_3
  end interface postpr_int_to_real

  interface postpr_read_restart
     module procedure &
          postpr_read_restart_RP_1  , postpr_read_restart_RP_2  , postpr_read_restart_RP_3,&
          postpr_read_restart_IP_1  , postpr_read_restart_IP_2  , &
          postpr_read_restart_R3P_1
  end interface postpr_read_restart

  interface postpr_write_restart
     module procedure &
          postpr_write_restart_RP_1  , postpr_write_restart_RP_2 , postpr_write_restart_RP_3,&
          postpr_write_restart_IP_1  , postpr_write_restart_IP_2 , &
          postpr_write_restart_R3P_1 
  end interface postpr_write_restart


  interface postpr_right_now
     module procedure &
          postpr_right_now_i1,&
          postpr_right_now_r1,&
          postpr_right_now_r2,&
          postpr_right_now_r0
  end interface postpr_right_now

  public :: postpr                        ! Postprocess a variable
  public :: postpr_postprocess            ! Postprocess a using mod_arrays format
  public :: postpr_read_restart           ! Read the restart of a variable
  public :: postpr_write_restart          ! Write the restart of a variable
  public :: postpr_right_now              ! Postprocess a variable on the fly
  public :: postpr_right_now_r0           ! Postprocess a variable on the fly
  public :: postpr_at_current_time_step   ! Postprocess at current step
  public :: postpr_initialization         ! Initialize postprocess variables

contains
 
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Initialize module variables
  !> @details Initialize variables and nullify arrays of this module
  !>
  !----------------------------------------------------------------------

  subroutine postpr_initialization()

    iiiii = 0_4
    rrrrr = 0.0_8
    wwwww = ''
    wwww8 = ''
    nullify(gescp)
    nullify(gevep)
    nullify(giscp)
    nullify(givep)
    nullify(bridge_ip_vector)
    nullify(bridge_ip_scalar)
    nullify(bridge_rp_vector)
    nullify(bridge_rp_scalar)

  end subroutine postpr_initialization

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of real scalar variables
  !> @details Bridge to the postprocess of a real variables
  !>
  !----------------------------------------------------------------------

  subroutine postpr_int_to_real_1(xx_int,xx_real)

    integer(ip), intent(in),     pointer :: xx_int(:)
    real(rp),    intent(inout),  pointer :: xx_real(:)
    integer(ip)                          :: ndim1
    integer(ip)                          :: idim1

    ndim1 = memory_size(xx_int)
    if( .not. associated(xx_real) ) then
       call memory_alloca(memke,'XX_REAL','mod_postpr',xx_real,ndim1)
    end if
    do idim1 = 1,ndim1
       xx_real(idim1) = real(xx_int(idim1),rp)
    end do
    
  end subroutine postpr_int_to_real_1
  
  subroutine postpr_int_to_real_2(xx_int,xx_real)

    integer(ip), intent(in),     pointer :: xx_int(:,:)
    real(rp),    intent(inout),  pointer :: xx_real(:,:)
    integer(ip)                          :: ndim1,ndim2
    integer(ip)                          :: idim1,idim2

    ndim1 = memory_size(xx_int,1_ip)
    ndim2 = memory_size(xx_int,2_ip)
    if( .not. associated(xx_real) ) then
       call memory_alloca(memke,'XX_REAL','mod_postpr',xx_real,ndim1,ndim2)
    end if
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          xx_real(idim1,idim2) = real(xx_int(idim1,idim2),rp)
       end do
    end do
    
  end subroutine postpr_int_to_real_2
  
  subroutine postpr_int_to_real_3(xx_int,xx_real)

    integer(ip), intent(in),     pointer :: xx_int(:,:,:)
    real(rp),    intent(inout),  pointer :: xx_real(:,:,:)
    integer(ip)                          :: ndim1,ndim2,ndim3
    integer(ip)                          :: idim1,idim2,idim3

    ndim1 = memory_size(xx_int,1_ip)
    ndim2 = memory_size(xx_int,2_ip)
    ndim3 = memory_size(xx_int,3_ip)
    if( .not. associated(xx_real) ) then
       call memory_alloca(memke,'XX_REAL','mod_postpr',xx_real,ndim1,ndim2,ndim3)
    end if
    do idim3 = 1,ndim3
       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             xx_real(idim1,idim2,idim3) = real(xx_int(idim1,idim2,idim3),rp)
          end do
       end do
    end do
    
  end subroutine postpr_int_to_real_3
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess "sobre la marcha"
  !> @details Postprocess a variable on the fly
  !>
  !----------------------------------------------------------------------

  subroutine postpr_right_now_i1(wopo1,wopo2,wopo3,iscalar)
    character(5),         intent(in) :: wopo1,wopo2,wopo3
    integer(ip), pointer, intent(in) :: iscalar(:)
    character(5)                     :: wopos_loc(3)
    integer(ip)                      :: ipoin,nbound

    if( INOTMASTER ) then
       if( wopo3 == 'NELEM') then
          nbound = nelem
       else
          nbound = npoin
       end if
       call postpr_memory(0_ip,nbound,0_ip)
       do ipoin = 1,nbound
          gescp(ipoin) = real(iscalar(ipoin),rp)
       end do
    end if
    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    call postpr(gescp,wopos_loc,0_ip,0.0_rp)
    if( INOTMASTER ) call postpr_memory(2_ip,nbound,0_ip)

  end subroutine postpr_right_now_i1

  subroutine postpr_right_now_r1(wopo1,wopo2,wopo3,rvector,pleng_opt)
    character(5),         intent(in)             :: wopo1,wopo2,wopo3
    real(rp),    pointer, intent(inout)          :: rvector(:)
    integer(ip),          intent(in),   optional :: pleng_opt
    character(5)                                 :: wopos_loc(3)

    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    call postpr(rvector,wopos_loc,0_ip,0.0_rp,pleng_opt=pleng_opt)

  end subroutine postpr_right_now_r1

  subroutine postpr_right_now_r0(wopo1,wopo2,wopo3,kpoin,rvector,TIME_STEP_OUTPUT,TIME_OUTPUT)
    integer(ip),  intent(in)             :: kpoin
    character(5), intent(in)             :: wopo1,wopo2,wopo3
    real(rp),     intent(inout)          :: rvector(*)
    real(rp),     intent(in),   optional :: TIME_OUTPUT
    integer(ip),  intent(in),   optional :: TIME_STEP_OUTPUT
    character(5)                         :: wopos_loc(3)
    integer(ip)                          :: ipoin
    real(rp)                             :: my_time
    integer(ip)                          :: my_time_step
    real(rp),     pointer                :: generic_scalar(:)

    nullify(generic_scalar)

    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    if( INOTMASTER ) then
       allocate(generic_scalar(npoin))
       do ipoin = 1,npoin
          generic_scalar(ipoin) = rvector(ipoin)
       end do
    end if
    if( present(TIME_OUTPUT) ) then
       my_time = TIME_OUTPUT
    else
       my_time = 0.0_rp
    end if
    if( present(TIME_STEP_OUTPUT) ) then
       my_time_step = TIME_STEP_OUTPUT
    else
       my_time_step = 0
    end if

    call postpr(generic_scalar,wopos_loc,my_time_step,my_time,TAG1=TIME_STEP_OUTPUT)

    if( INOTMASTER ) deallocate(generic_scalar)

  end subroutine postpr_right_now_r0

  subroutine postpr_right_now_r2(wopo1,wopo2,wopo3,rvector)
    character(5),         intent(in) :: wopo1,wopo2,wopo3
    real(rp),    pointer, intent(in) :: rvector(:,:)
    character(5)                     :: wopos_loc(3)
    integer(ip)                      :: ipoin,nbound
    real(rp),    pointer             :: my_vect(:,:)

    nullify(my_vect)
    if( INOTMASTER ) then
       if( wopo3 == 'NELEM') then
          nbound = nelem
       else
          nbound = npoin
       end if
       allocate(my_vect(ndime,nbound))
       do ipoin = 1,nbound
          my_vect(1:ndime,ipoin) = rvector(1:ndime,ipoin)
       end do
    end if
    wopos_loc(1) = wopo1
    wopos_loc(2) = wopo2
    wopos_loc(3) = wopo3
    call postpr(my_vect,wopos_loc,0_ip,0.0_rp)
    if( INOTMASTER ) deallocate(my_vect)

  end subroutine postpr_right_now_r2

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess current time step
  !> @details Postprocess current time step
  !>
  !----------------------------------------------------------------------

  function postpr_at_current_time_step()
    integer(ip) :: postpr_at_current_time_step
    integer(ip) :: imodu,ivari,itime,iok
    logical(lg) :: if_already

    postpr_at_current_time_step = 0

    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          do ivari = 1,nvarp

             if_already = .false.

             if( ittyp == ITASK_ENDRUN ) then
                iok = 1
             else
                iok = 0
             end if

             if( &
                  ittim >= momod(imodu) % postp(1) % npp_inits .and. &
                  iok == 0 .and. &
                  momod(imodu) % postp(1) % npp_stepi(ivari,1) > 0 ) then
                if( mod(ittim, momod(imodu) % postp(1) % npp_stepi(ivari,1)) == 0 ) then
                   postpr_at_current_time_step = postpr_at_current_time_step + 1
                   if_already = .true.
                end if
             end if
             !
             ! At a given time
             !
             if( iok == 0 .and. ( .not. if_already ) ) then
                do itime = 1,nvart
                   if(   abs( momod(imodu) % postp(1) % pos_times(itime,ivari,1)-cutim) < (0.5_rp*dtime) .and. &
                        &     momod(imodu) % postp(1) % pos_times(itime,ivari,1)        > 0.0_rp) then
                      postpr_at_current_time_step = postpr_at_current_time_step + 1
                      if_already = .true.
                   end if
                end do
             end if
             !
             ! At a given time period
             !
             if( cutim >= momod(imodu) % postp(1) % pos_tinit .and. iok == 0 &
                  .and. ( .not. if_already ) ) then
                if(    abs(momod(imodu) % postp(1) % pos_times(1,ivari,1)-cutim) < (0.6_rp*dtime).and.&
                     &     momod(imodu) % postp(1) % pos_perio(ivari,1)          > 0.0_rp) then
                   postpr_at_current_time_step = postpr_at_current_time_step + 1
                   if_already = .true.
                end if
             end if

          end do
       end if
    end do
  end function postpr_at_current_time_step
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of real vector variables
  !> @details Bridge to the postprocess of a real variables
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine postpr_postprocess_IP_1(bridge,ivari,itste,ttime,mesh)

    integer(ip),                     pointer, intent(inout) :: bridge(:)
    integer(ip),                              intent(in)    :: ivari
    integer(ip) ,                             intent(in)    :: itste
    real(rp)    ,                             intent(in)    :: ttime
    type(mesh_type_basic),  optional,         intent(in)    :: mesh
    real(rp),                        pointer                :: bridge_loc(:)

    nullify(bridge_loc)
    call postpr_int_to_real(bridge,bridge_loc)
    call postpr_postprocess(bridge_loc,ivari,itste,ttime,mesh)
    call memory_deallo(memke,'XX_REAL','mod_postpr',bridge_loc)

  end subroutine postpr_postprocess_IP_1
  
  subroutine postpr_postprocess_IP_2(bridge,ivari,itste,ttime,mesh)

    integer(ip),                     pointer, intent(inout) :: bridge(:,:)
    integer(ip),                              intent(in)    :: ivari
    integer(ip) ,                             intent(in)    :: itste
    real(rp)    ,                             intent(in)    :: ttime
    type(mesh_type_basic),  optional,         intent(in)    :: mesh
    real(rp),                        pointer                :: bridge_loc(:,:)

    nullify(bridge_loc)
    call postpr_int_to_real(bridge,bridge_loc)
    call postpr_postprocess(bridge_loc,ivari,itste,ttime,mesh)
    call memory_deallo(memke,'XX_REAL','mod_postpr',bridge_loc)

  end subroutine postpr_postprocess_IP_2
  
  subroutine postpr_postprocess_IP_3(bridge,ivari,itste,ttime,mesh)

    integer(ip),                     pointer, intent(inout) :: bridge(:,:,:)
    integer(ip),                              intent(in)    :: ivari
    integer(ip) ,                             intent(in)    :: itste
    real(rp)    ,                             intent(in)    :: ttime
    type(mesh_type_basic),  optional,         intent(in)    :: mesh
    real(rp),                        pointer                :: bridge_loc(:,:,:)

    nullify(bridge_loc)
    call postpr_int_to_real(bridge,bridge_loc)
    call postpr_postprocess(bridge_loc,ivari,itste,ttime,mesh)
    call memory_deallo(memke,'XX_REAL','mod_postpr',bridge_loc)

  end subroutine postpr_postprocess_IP_3
  
  subroutine postpr_postprocess_RP_1(bridge,ivari,itste,ttime,mesh,inte)

    real(rp),                        pointer, intent(inout) :: bridge(:)
    integer(ip),                              intent(in)    :: ivari
    integer(ip) ,                             intent(in)    :: itste
    real(rp)    ,                             intent(in)    :: ttime
    type(mesh_type_basic),  optional,         intent(in)    :: mesh
    type(typ_interp),       optional,         intent(in)    :: inte
    real(rp),     pointer                                   :: vector_loc(:,:)
    real(rp),     pointer                                   :: scalar_loc(:)

    nullify(scalar_loc,vector_loc)
    call postpr_options(ivari)
    if( wopos_loc(2) == 'VECTO' ) then
       call runend('postpr_RP_1: CANNOT POSTPROCESS A VECTOR AS A SCALAR')
    else if( wopos_loc(2) == 'SCALA' ) then
       if( present(inte) .and. present(mesh) ) then
          nullify(xx_sca)
          call interp_fe(bridge,inte,xx_sca,mesh % npoin)
          call postpr_real_scalar(xx_sca,wopos_loc,itste,ttime,MESH=mesh)
          call interp_fe_deallocate(xx_sca)
       else                
          call postpr_real_scalar(bridge,wopos_loc,itste,ttime,MESH=mesh)
       end if
    else
       call runend('postpr_RP_1: NO NOT KNOW WHAT TO DO')
    end if

  end subroutine postpr_postprocess_RP_1
  
  subroutine postpr_postprocess_RP_2(bridge,ivari,itste,ttime,mesh,inte)

    real(rp),                        pointer, intent(inout) :: bridge(:,:)
    integer(ip),                              intent(in)    :: ivari
    integer(ip) ,                             intent(in)    :: itste
    real(rp)    ,                             intent(in)    :: ttime
    type(mesh_type_basic),  optional,         intent(in)    :: mesh
    type(typ_interp),       optional,         intent(in)    :: inte
    integer(ip)                                             :: ienti,icomp
    integer(ip)                                             :: pdime,pcomp,kcomp   
    real(rp),     pointer                                   :: vector_loc(:,:)
    real(rp),     pointer                                   :: scalar_loc(:)
    character(5)                                            :: wopoc
    
    nullify(scalar_loc,vector_loc)
    call postpr_options(ivari,MESH=mesh)
    if( wopos_loc(2) == 'SCALA' ) then
       call memory_alloca(memke,'SCALAR_LOC','mod_postpr',scalar_loc,nenti)
    end if

    if( comp_posit_loc == 0 .and. wopos_loc(2) == 'VECTO' ) then
       !
       ! This is a vector without components
       !
       if( enti_posit_loc == 2 ) then
          if( postp(1) % dime_num(ivari) == 0 ) then
             pdime = ndime
          else
             pdime = postp(1) % dime_num(ivari)
          end if
          if( present(inte) .and. present(mesh) ) then
             nullify(xx_vec)
             call interp_fe(bridge,inte,xx_vec,mesh % npoin)
             call postpr_real_vector(bridge,wopos_loc,itste,ttime,pdime,MESH=mesh)
             call interp_fe_deallocate(xx_vec)
          else
             call postpr_real_vector(bridge,wopos_loc,itste,ttime,pdime,MESH=mesh)
          end if
       else
          call runend('postpr_real_RP_3: DO NOT KNOW WHAT TO DO FOR '//wopos_loc(1))
       end if
       
    else if( comp_posit_loc == 0 .and. wopos_loc(2) == 'SCALA' ) then
       !
       ! This is a scalar without components
       !
       if( enti_posit_loc == 1 ) then
          do ienti = 1,nenti
             scalar_loc(ienti) = bridge(ienti,1)             
          end do
       else if( enti_posit_loc == 2 ) then
          do ienti = 1,nenti
             scalar_loc(ienti) = bridge(1,ienti)             
          end do          
       else
          call runend('postpr_real_RP_3: DO NOT KNOW WHAT TO DO FOR '//wopos_loc(1))
       end if
       if( present(inte) .and. present(mesh) ) then
          nullify(xx_sca)
          call interp_fe(scalar_loc,inte,xx_sca,mesh % npoin)
          call postpr_real_scalar(xx_sca,wopos_loc,itste,ttime,MESH=mesh)
          call interp_fe_deallocate(xx_sca)
       else          
          call postpr_real_scalar(scalar_loc,wopos_loc,itste,ttime,MESH=mesh)
       end if
       
    else if( comp_posit_loc > 0 .and. wopos_loc(2) == 'SCALA' ) then
       !
       ! This is a scalar with component
       !
       if( postp(1) % lcomp(1,nvarp,id_loc) == -1 ) then
          pcomp = postp(1) % comp_num(ivari)
          if( pcomp == 0 ) then
             !
             ! Composant nuber could not be guessed when allocating...
             ! this is probably a secondary variable
             !
             pcomp = memory_size(bridge,comp_posit_loc)
             call PAR_MAX(pcomp)
             postp(1) % comp_num(ivari) = pcomp
          end if
       else
          pcomp = size(postp(1) % lcomp,1_ip)
       end if
       
       do kcomp = 1,pcomp

          if( postp(1) % lcomp(1,ivari,id_loc) == -1 ) then
             icomp = kcomp
          else
             icomp = postp(1) % lcomp(kcomp,ivari,id_loc) 
          end if

          if( icomp > 0 ) then
             wopoc = wopos_loc(1)
             call postpr_components(icomp,wopoc,postp(1) % comp_num(ivari))
             wopos_loc(1) = wopoc
             if( comp_posit_loc == 1 ) then
                do ienti = 1,nenti
                   scalar_loc(ienti) = bridge(icomp,ienti)
                end do
             else if( comp_posit_loc == 2 ) then
                do ienti = 1,nenti
                   scalar_loc(ienti) = bridge(ienti,icomp)
                end do
             end if
             if( present(inte) .and. present(mesh) ) then
                nullify(xx_sca)
                call interp_fe(scalar_loc,inte,xx_sca,mesh % npoin)
                call postpr_real_scalar(xx_sca,wopos_loc,itste,ttime,MESH=mesh)
                call interp_fe_deallocate(xx_sca)
             else          
                call postpr_real_scalar(scalar_loc,wopos_loc,itste,ttime,MESH=mesh)
             end if
          end if
          
       end do
    end if
    
    if( wopos_loc(2) == 'SCALA' ) then
       call memory_deallo(memke,'SCALAR_LOC','mod_postpr',scalar_loc)
    end if
    
  end subroutine postpr_postprocess_RP_2
  
  subroutine postpr_postprocess_RP_3(bridge,ivari,itste,ttime,mesh,inte)

    real(rp),                        pointer, intent(inout) :: bridge(:,:,:)
    integer(ip),                              intent(in)    :: ivari
    integer(ip) ,                             intent(in)    :: itste
    real(rp)    ,                             intent(in)    :: ttime
    type(mesh_type_basic),  optional,         intent(in)    :: mesh
    type(typ_interp),       optional,         intent(in)    :: inte
    integer(ip)                                             :: ienti,idime,icomp
    integer(ip)                                             :: pdime,pcomp,kcomp    
    real(rp),     pointer                                   :: vector_loc(:,:)
    real(rp),     pointer                                   :: scalar_loc(:)
    character(5)                                            :: wopoc
    
    nullify(scalar_loc,vector_loc)
    call postpr_options(ivari,MESH=mesh)
    
    if( wopos_loc(2) == 'SCALA' ) then
       call memory_alloca(memke,'SCALAR_LOC','mod_postpr',scalar_loc,nenti)
    else if( wopos_loc(2) == 'VECTO' ) then
       call memory_alloca(memke,'VECTOR_LOC','mod_postpr',vector_loc,pdime,nenti)
    end if

    if( comp_posit_loc > 0 .and. wopos_loc(2) == 'VECTO' ) then
       !
       ! This is a vector with components
       !
       call runend('postpr_real_RP_3: VECTOR WITH COMPONENTS NOT CODED')

    else if( comp_posit_loc == 0 .and. wopos_loc(2) == 'VECTO' ) then
       !
       ! This is a vector without components
       !
       if( enti_posit_loc == 2 ) then
          if( postp(1) % dime_num(ivari) == 0 ) then
             pdime = ndime
          else
             pdime = postp(1) % dime_num(ivari)
          end if
          do ienti = 1,nenti
             do idime = 1,pdime
                vector_loc(idime,ienti) = bridge(idime,ienti,1)
             end do
          end do
          if( present(inte) .and. present(mesh) ) then
             nullify(xx_vec)
             call interp_fe(vector_loc,inte,xx_vec,mesh % npoin)
             call postpr_real_vector(xx_vec,wopos_loc,itste,ttime,pdime,MESH=mesh)
             call interp_fe_deallocate(xx_vec)
          else
             call postpr_real_vector(vector_loc,wopos_loc,itste,ttime,pdime,MESH=mesh)
          end if
       else
          call runend('postpr_real_RP_3: DO NOT KNOW WHAT TO DO FOR '//wopos_loc(1))
       end if
       
    else if( comp_posit_loc == 0 .and. wopos_loc(2) == 'SCALA' ) then
       !
       ! This is a scalar without component
       !
       if( enti_posit_loc == 1 ) then
          
          do ienti = 1,nenti
             scalar_loc(ienti) = bridge(ienti,1,1)
          end do
          
       else if( enti_posit_loc == 2 ) then
          
          do ienti = 1,nenti
             scalar_loc(ienti) = bridge(1,ienti,1)
          end do

       else 
          
          do ienti = 1,nenti
             scalar_loc(ienti) = bridge(1,1,ienti)
          end do

       end if
       if( present(inte) .and. present(mesh) ) then
          nullify(xx_sca)
          call interp_fe(scalar_loc,inte,xx_sca,mesh % npoin)
          call postpr_real_scalar(xx_sca,wopos_loc,itste,ttime,MESH=mesh)
          call interp_fe_deallocate(xx_sca)
       else          
          call postpr_real_scalar(scalar_loc,wopos_loc,itste,ttime,MESH=mesh)
       end if
       
    else if( comp_posit_loc > 0 .and. wopos_loc(2) == 'SCALA' ) then
       !
       ! This is a scalar with component
       !       
       if( postp(1) % lcomp(1,ivari,id_loc) == -1 ) then
          pcomp = postp(1) % comp_num(ivari)
       else
          pcomp = size(postp(1) % lcomp,1_ip)
       end if

       do kcomp = 1,pcomp

          if( postp(1) % lcomp(1,ivari,id_loc) == -1 ) then
             icomp = kcomp
          else
             icomp = postp(1) % lcomp(kcomp,ivari,id_loc) 
          end if
          
          if( icomp > 0 ) then
             wopoc = wopos_loc(1)
             call postpr_components(icomp,wopoc,postp(1) % comp_num(ivari))
             wopos_loc(1) = wopoc
             if( comp_posit_loc == 1 .and. enti_posit_loc == 2 ) then
                do ienti = 1,nenti
                   scalar_loc(ienti) = bridge(icomp,ienti,1)
                end do
             else if( comp_posit_loc == 2 .and. enti_posit_loc == 1 ) then
                do ienti = 1,nenti
                   scalar_loc(ienti) = bridge(ienti,icomp,1)
                end do
             else
                call runend('postpr_real_RP_3: DO NOT KNOW WHAT TO DO FOR '//wopos_loc(1))
             end if
             
          end if
          if( present(inte) .and. present(mesh) ) then
             nullify(xx_sca)
             call interp_fe(scalar_loc,inte,xx_sca,mesh % npoin)
             call postpr_real_scalar(xx_sca,wopos_loc,itste,ttime,MESH=mesh)
             call interp_fe_deallocate(xx_sca)
          else          
             call postpr_real_scalar(scalar_loc,wopos_loc,itste,ttime,MESH=mesh)
          end if
       end do
    end if

    if( wopos_loc(2) == 'SCALA' ) then
       call memory_deallo(memke,'SCALAR_LOC','mod_postpr',scalar_loc)
    else if( wopos_loc(2) == 'VECTO' ) then
       call memory_deallo(memke,'VECTOR_LOC','mod_postpr',vector_loc)
    end if
    
  end subroutine postpr_postprocess_RP_3
  
  subroutine postpr_postprocess_R3P_1(bridge,ivari,itste,ttime)

    type(r3p),    intent(inout), pointer :: bridge(:)
    integer(ip),  intent(in)             :: ivari
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime

    call postpr_options(ivari)
    call postpr(bridge,wopos_loc,itste,ttime)
    
  end subroutine postpr_postprocess_R3P_1
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of real scalar variables
  !> @details Bridge to the postprocess of a real variables
  !>
  !----------------------------------------------------------------------

  subroutine postpr_real_scalar(bridge,wopos,itste,ttime,worig,pleng_opt,TAG1,TAG2,mesh)

    real(rp),                         pointer, intent(inout) :: bridge(:)
    character(5),                              intent(in)    :: wopos(*)
    integer(ip),                               intent(in)    :: itste
    real(rp),                                  intent(in)    :: ttime
    character(5),           optional,          intent(in)    :: worig
    integer(ip),            optional,          intent(in)    :: pleng_opt
    integer(ip),            optional,          intent(in)    :: TAG1
    integer(ip),            optional,          intent(in)    :: TAG2
    type(mesh_type_basic),  optional,          intent(in)    :: mesh
    integer(ip)                                              :: pdime

    pdime     = 1
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if

    if( present(mesh) ) then
       !
       ! Postprocess on a mesh
       !
       if (postpr_is_mpio(wopos,itste,ttime,pdime,0_ip,TAG1,TAG2,.false.)) then
          call postpr_mesh_mpio(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       else
          call postpr_mesh_alya(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       end if
       
    else if (postpr_is_mpio(wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)) then
       !
       ! MPIO
       !
       gescar_mpio => bridge
       call posmpio_real_v()

    else
       !
       ! Alya
       !
       if( kfl_reawr /= 1 ) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)

       !if ( kfl_reawr /= 1 ) then
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2)
       !else
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2,READ_MODE=.true.)
       !end if
       if( IMASTER ) then
          call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,mesh)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_rp(bridge,    wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,mesh)
          else
             call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,mesh)
          end if
       end if
       call end_timer()

    end if

  end subroutine postpr_real_scalar
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of real vector variables
  !> @details Bridge to the postprocess of a real variables
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine postpr_real_matrix(bridge,wopos,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2)

    implicit none
    character(5), intent(in)             :: wopos(*)
    real(rp),     intent(inout), pointer :: bridge(:,:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  intent(in), optional   :: kdime
    character(5), intent(in), optional   :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip)                          :: idime,nenti,ndim1,ndim2
    integer(ip)                          :: idim1,idim2,ienti,pdime
    !
    ! The master needs the dimension
    !
    ndim1 = memory_size(bridge,1_ip)
    ndim2 = memory_size(bridge,2_ip) 
    nenti = memory_size(bridge,3_ip)
    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = ndim1*ndim2
       call PAR_MAX(pdime)
    end if
    call memory_alloca(memke,wopos(1),'mod_postpr',bridge_rp_vector,pdime,nenti)

    if( kfl_reawr /= 1 ) then
       do ienti  = 1,nenti
          idime = 0
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                idime = idime + 1
                bridge_rp_vector(idime,ienti) = bridge(idim1,idim2,ienti) 
             end do
          end do
       end do
    end if

    call postpr_real_vector(bridge_rp_vector,wopos,itste,ttime,pdime,worig,pleng_opt,TAG1,TAG2)

    if( kfl_reawr == 1 ) then
       do ienti  = 1,nenti
          idime = 0
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                idime = idime + 1
                bridge(idim1,idim2,ienti) = bridge_rp_vector(idime,ienti)
             end do
          end do
       end do
    end if

    call memory_deallo(memke,wopos(1),'mod_postpr',bridge_rp_vector)

  end subroutine postpr_real_matrix

  subroutine postpr_real_vector(bridge,wopos,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,mesh)

    real(rp),                         pointer, intent(inout) :: bridge(:,:)
    character(5),                              intent(in)    :: wopos(*)
    integer(ip),                               intent(in)    :: itste
    real(rp),                                  intent(in)    :: ttime
    integer(ip),            optional,          intent(in)    :: kdime
    character(5),           optional,          intent(in)    :: worig
    integer(ip),            optional,          intent(in)    :: pleng_opt
    integer(ip),            optional,          intent(in)    :: TAG1
    integer(ip),            optional,          intent(in)    :: TAG2
    type(mesh_type_basic),  optional,          intent(in)    :: mesh
    integer(ip)                                              :: pdime,dummi
    
    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = ndime
    end if
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if

    if( present(mesh) ) then
       !
       ! Postprocess on a mesh
       !
       dummi = 0_ip
       if (postpr_is_mpio(wopos,itste,ttime,pdime,dummi,TAG1,TAG2,.false.)) then
          call postpr_mesh_mpio(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       else
          call postpr_mesh_alya(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       end if
       
    else if (postpr_is_mpio(wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2)) then
       !
       ! MPIO
       !
       gevecr_mpio    => bridge
       call posmpio_real_m()

    else
       !
       ! Alya
       !
       if( kfl_reawr /=1 ) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)
       !if ( kfl_reawr /= 1 ) then
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2)
       !else
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2,READ_MODE=.true.)
       !end if
       if( IMASTER ) then
          call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,mesh)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_rp(bridge,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,mesh)
          else
             call postpr_postprocess_rp(bridge_rp,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,mesh)
          end if
       end if
       call end_timer()

    end if

  end subroutine postpr_real_vector

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of integer scalar variables
  !> @details Bridge to the postprocess of a integer variables
  !>
  !----------------------------------------------------------------------

  subroutine postpr_int_scalar(bridge,wopos,itste,ttime,worig,TAG1,TAG2,mesh)

    integer(ip),                      pointer, intent(inout) :: bridge(:)
    character(5),                              intent(in)    :: wopos(*)
    integer(ip),                               intent(in)    :: itste
    real(rp),                                  intent(in)    :: ttime
    character(5),           optional,          intent(in)    :: worig
    integer(ip),            optional,          intent(in)    :: TAG1
    integer(ip),            optional,          intent(in)    :: TAG2
    type(mesh_type_basic),  optional,          intent(in)    :: mesh
    integer(ip)                                              :: pdime

    pdime = 1
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if
    
    if( present(mesh) ) then
       !
       ! Postprocess on a mesh
       !
       if (postpr_is_mpio(wopos,itste,ttime,pdime,0_ip,TAG1,TAG2,.false.)) then
          call postpr_mesh_mpio(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       else
          call postpr_mesh_alya(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       end if
       
    else if (postpr_is_mpio(wopos,itste,ttime,pdime,0_ip,TAG1,TAG2)) then
       !
       ! MPIO
       !
       gescai_mpio    => bridge
       call posmpio_int_v()

    else
       !
       ! Alya
       !
       if( kfl_reawr /= 1 ) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)

       !if ( kfl_reawr /= 1 ) then
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2)
       !else
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2,READ_MODE=.true.)
       !end if

       if( IMASTER ) then
          call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2,mesh)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_ip(bridge,wopos,itste,ttime,pdime,TAG1,TAG2,mesh)
          else
             call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2,mesh)
          end if
       end if
       call end_timer()
    end if

#ifdef EVENT
    call mpitrace_user_function(0)
#endif

  end subroutine postpr_int_scalar

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of integer vector variables
  !> @details Bridge to the postprocess of a integer variables
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine postpr_int_vector(bridge,wopos,itste,ttime,kdime,worig,TAG1,TAG2,mesh)

    integer(ip),                      pointer, intent(inout) :: bridge(:,:)
    character(5),                              intent(in)    :: wopos(*)
    integer(ip),                               intent(in)    :: itste
    real(rp),                                  intent(in)    :: ttime
    integer(ip),            optional,          intent(in)    :: kdime
    character(5),           optional,          intent(in)    :: worig
    integer(ip),            optional,          intent(in)    :: TAG1
    integer(ip),            optional,          intent(in)    :: TAG2
    type(mesh_type_basic),  optional,          intent(in)    :: mesh
    integer(ip)                                              :: pdime

    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = ndime
    end if
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if

    if( present(mesh) ) then
       !
       ! Postprocess on a mesh
       !
       if (postpr_is_mpio(wopos,itste,ttime,pdime,0_ip,TAG1,TAG2,.false.)) then
          call postpr_mesh_mpio(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       else
          call postpr_mesh_alya(mesh,bridge,wopos,itste,ttime,TAG1,TAG2)
       end if
       
    else if (postpr_is_mpio(wopos,itste,ttime,pdime,0_ip,TAG1,TAG2)) then
       !
       ! MPIO
       !
       geveci_mpio => bridge
       call posmpio_int_m()
       
    else
       !
       ! Alya
       !
       if( kfl_reawr /= 1 ) then
          bench_rw='WRITE'
       else
          bench_rw='READ '
       end if
       call start_timer(wopos(1), bench_rw, 'ALYA', 1_ip, 1_ip)
       !if ( kfl_reawr /= 1 ) then
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2)
       !else
       !   call postpr_mesh_alya(meshe(ndivi),bridge,wopos,itste,ttime,TAG1,TAG2,READ_MODE=.true.)
       !end if
       if( IMASTER ) then
          call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2,mesh)
       else
          if( associated(bridge) ) then
             call postpr_postprocess_ip(bridge,wopos,itste,ttime,pdime,TAG1,TAG2,mesh)
          else
             call postpr_postprocess_ip(bridge_ip,wopos,itste,ttime,pdime,TAG1,TAG2,mesh)
          end if
       end if
       call end_timer()
    end if

  end subroutine postpr_int_vector

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of r3p type
  !> @details Bridge to the postprocess of a r3p types
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine posr3p(bridge,wopos,itste,ttime,kdime,TAG1,TAG2)

    implicit none
    character(5), intent(in)              :: wopos(*)
    type(r3p),    intent(inout), pointer  :: bridge(:)
    integer(ip),  intent(in)              :: itste
    real(rp),     intent(in)              :: ttime
    integer(ip),  intent(in),    optional :: kdime
    integer(ip),  intent(in),    optional :: TAG1
    integer(ip),  intent(in),    optional :: TAG2
    integer(ip)                           :: ileng,pleng
    integer(ip)                           :: idim1,idim2,ptota,ielty,ii
    integer(ip)                           :: pdim1,pdim2,itota,pdime

    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = 1
    end if
       !
       ! Type of postprocess
       !
       wwwww(1) = 'ALYA '
       wwwww(2) = 'V0003'
       wwwww(3) = wopos(1)
       wwwww(4) = wopos(2)
       wwwww(5) = 'R3P  '
       wwwww(6) = 'REAL '
       wwwww(7) = '8BYTE'
       if( wopos(2) == 'R3PVE' ) then
          iiiii(1) = int(ndime,4)
       else
          iiiii(1) = 1_4
       end if
       iiiii(2) = 0_4
       iiiii(4) = int(itste,4)
       rrrrr(1) = ttime
       wopos_pos(1) = wopos(1)
       parii    = NELEM_TYPE
       pleng    = nelem
       if( IMASTER ) pleng = nelem_total

       iiiii(2) = int(pleng,4)
       !
       ! Output solution
       !
       call postpr_open_file(TAG1,TAG2)
       !call opfpos(1_ip)
       if( kfl_reawr < 0 ) then
          kfl_reawr = abs(kfl_reawr)
          return
       end if

       call postpr_header()

       if( ISEQUEN ) then

          if( kfl_reawr /= 1 ) then
             ptota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                ptota = ptota + pdim1 * pdim2
             end do
             write(lun_postp) ( ngaus(ielty),ielty=iesta_dom,iesto_dom)
             write(lun_postp) pleng
             write(lun_postp) ptota
             write(lun_postp) (((bridge(ileng)%a(idim1,idim2,pdime),idim1=1,size(bridge(ileng)%a,1,kind=ip)),&
                  idim2=1,size(bridge(ileng)%a,2,kind=ip)),ileng=1,pleng)
          else
             read(lun_postp)  ( ii,ielty=iesta_dom,iesto_dom)
             read(lun_postp)  pleng
             read(lun_postp)  ptota
             read(lun_postp)  (((bridge(ileng)%a(idim1,idim2,pdime),idim1=1,size(bridge(ileng)%a,1,kind=ip)),&
                  idim2=1,size(bridge(ileng)%a,2,kind=ip)),ileng=1,pleng)
          end if

       else if( IMASTER ) then

          if( kfl_reawr /= 1 ) then
             write(lun_postp) ( ngaus(ielty),ielty=iesta_dom,iesto_dom)
          else if( kfl_reawr == 1 ) then
             read(lun_postp) ( ii,ielty=iesta_dom,iesto_dom)
          end if

          do kfl_desti_par = 1,npart

             if( kfl_reawr /= 1 ) then
                !
                ! - Master writes postpro without filter
                ! - Master writes restart
                !
                pleng = nelem_par(kfl_desti_par)
                call pararr('RCV',0_ip,1_ip,ptota)
                write(lun_postp) pleng
                write(lun_postp) ptota
                call postpr_memory(0_ip,ptota,0_ip)
                call pararr('RCV',0_ip,ptota,gescp)
                write(lun_postp) ( gescp(ileng), ileng = 1,ptota )
                call postpr_memory(2_ip,ptota,0_ip)

             else if( kfl_reawr == 1 ) then
                !
                ! - Master reads restart and sends to slave
                !
                pleng = nelem_par(kfl_desti_par)
                read(lun_postp) pleng
                read(lun_postp) ptota
                call postpr_memory(0_ip,ptota,0_ip)
                read(lun_postp) ( gescp(ileng), ileng = 1,ptota )
                call pararr('SND',0_ip,ptota,gescp)
                call postpr_memory(2_ip,ptota,0_ip)

             end if

          end do

       else if( ISLAVE ) then

          kfl_desti_par = 0

          if( kfl_reawr /= 1 ) then
             !
             ! Slaves without filter
             !
             ptota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                ptota = ptota + pdim1 * pdim2
             end do
             call parari('SND',0_ip,1_ip,ptota)
             call postpr_memory(0_ip,ptota,0_ip)
             itota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                do idim2=1,pdim2
                   do idim1=1,pdim1
                      itota = itota + 1
                      gescp(itota) = bridge(ileng)%a(idim1,idim2,pdime)
                   end do
                end do
             end do
             call pararr('SND',0_ip,ptota,gescp)
             call postpr_memory(2_ip,ptota,0_ip)

          else if( kfl_reawr == 1 ) then
             !
             ! Slaves receive
             !
             ptota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                ptota = ptota + pdim1 * pdim2
             end do
             call postpr_memory(0_ip,ptota,0_ip)
             call pararr('RCV',0_ip,ptota,gescp)
             itota = 0
             do ileng = 1,nelem
                pdim1 = size(bridge(ileng)%a,1,kind=ip)
                pdim2 = size(bridge(ileng)%a,2,kind=ip)
                do idim2=1,pdim2
                   do idim1=1,pdim1
                      itota = itota + 1
                      bridge(ileng)%a(idim1,idim2,pdime) = gescp(itota)
                   end do
                end do
             end do
             call postpr_memory(2_ip,ptota,0_ip)

          end if

       end if

       call postpr_close_file()

#ifdef EVENT
    call mpitrace_user_function(0)
#endif

  end subroutine posr3p

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess of complex vector variables
  !> @details Bridge to the postprocess of a complex variables
  !>          When kdime is present, kdime rules over ndime.
  !>
  !----------------------------------------------------------------------

  subroutine posvex(bridge,wopos,itste,ttime,kdime,worig)

    implicit none
    character(5), intent(in)             :: wopos(*)
    complex(rp),  intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip)                          :: pdime

    if( present(kdime) ) then
       pdime = kdime
    else
       pdime = ndime
    end if
    kfl_origi = 1
    if( present(worig) ) then
       if( worig == 'ORIGI' ) kfl_origi = 0
    end if

  end subroutine posvex

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess a real array
  !> @details Postprocess a real array
  !>
  !----------------------------------------------------------------------

  subroutine postpr_postprocess_rp(bridge,wopos,itste,ttime,pdime,pleng_opt,TAG1,TAG2,mesh)

    character(5),                      intent(in)    :: wopos(*)
    real(rp),                          intent(inout) :: bridge(*)
    integer(ip),                       intent(in)    :: itste
    real(rp),                          intent(in)    :: ttime
    integer(ip),                       intent(in)    :: pdime
    integer(ip),            optional,  intent(in)    :: pleng_opt
    integer(ip),            optional,  intent(in)    :: TAG1
    integer(ip),            optional,  intent(in)    :: TAG2
    type(mesh_type_basic),  optional,  intent(in)    :: mesh
    integer(ip)                                      :: ileng,pleng,dummi(2)
    integer(ip)                                      :: imesh,pleng_total
    integer(ip)                                      :: pleng_com(2),tag1_opt,tag2_opt

    !
    ! Options
    !
    tag1_opt = 0
    tag2_opt = 0
    if( present(TAG1) ) tag1_opt = TAG1
    if( present(TAG2) ) tag2_opt = TAG2
    !
    ! Type of postprocess
    !
    wwwww(1) = 'ALYA '
    wwwww(2) = 'V0003'
    wwwww(3) = wopos(1)
    wwwww(4) = wopos(2)
    wwwww(5) = wopos(3)
    iiiii(5) = int(tag1_opt,4)
    iiiii(6) = int(tag2_opt,4)
    wwwww(6) = 'REAL '
    wwwww(7) = '8BYTE'
    iiiii(1) = int(pdime,4)
    iiiii(2) = 0_4
    iiiii(4) = int(itste,4)
    rrrrr(1) = ttime

    wopos_pos(1) = wopos(1)

    if( kfl_reawr == 0 ) then
       imesh = kfl_posdi
    else
       imesh = ndivi
    end if
    !
    ! Map the result onto original mesh:
    ! - Only for NPOIN type results
    !
    kfl_permu = 0
    kfl_perme = 0
    kfl_permb = 0

    if( kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos(3) == 'NPOIN' ) kfl_permu = 1
    if( kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos(3) == 'NELEM' ) kfl_perme = 1
    if( kfl_origi == 1 .and. ( kfl_posdi /= ndivi ) .and. wopos(3) == 'NBOUN' ) kfl_permb = 1

    if( wopos(3) == 'NPOIN' ) then
       !
       ! NPOIN variable
       !
       parii = NPOIN_TYPE
       if( associated(meshe) .and. size(meshe,kind=ip) >= imesh ) then
          pleng = meshe(imesh) % npoin
          if( IMASTER ) pleng = meshe(imesh) % npoin_total
       else
          pleng = npoin
          if( IMASTER ) pleng = npoin_total
       end if
       
    else if( wopos(3) == 'NELEM' ) then
       !
       ! NELEM variable
       !
       parii = NELEM_TYPE
       pleng = meshe(imesh) % nelem
       if( IMASTER ) pleng = meshe(imesh) % nelem_total
       
    else if( wopos(3) == 'NBOUN' ) then
       !
       ! NBOUN variable
       !
       parii = NBOUN_TYPE
       pleng = meshe(imesh) % nboun
       if( IMASTER ) pleng = meshe(imesh) % nboun_total

    else if( wopos(3) == 'WHATE' ) then

       parii = NBOUN_TYPE
       if( present(pleng_opt) ) then
          pleng = pleng_opt
       else
          call runend('MOD_POSTPR - POSTPR_POSTPROCESS_RP: A SIZE SHOULD BE GIVEN')
       end if
       if( IMASTER ) then
          pleng_total = 0
       else
          pleng_total = pleng
       end if
       call PAR_SUM(pleng_total)
       if( IMASTER ) pleng = pleng_total

    else
       call runend('MOD_POSTPR - POSTPR_POSTPROCESS_RP: UNDEFINED POSTPROCESS')
    end if
    iiiii(2) = int(pleng,4)
    !
    ! Output solution
    !
    call postpr_open_file(TAG1,TAG2)

    if( kfl_reawr < 0 ) then
       kfl_reawr = abs(kfl_reawr)
       return
    end if
    !
    ! Write header
    !
    call postpr_header()
    
    if( ISEQUEN ) then

       !-----------------------------------------------------------------
       !
       ! SEQUENTIAL
       !
       !-----------------------------------------------------------------

       if( kfl_reawr == 0 .or. kfl_reawr == 2 ) then
          !
          ! Postprocess and preliminary 
          !
          write(lun_postp) pleng
          write(lun_postp) ( bridge(ileng), ileng = 1,pdime*pleng )

       else if( kfl_reawr == 1 ) then
          !
          ! Restart
          !
          read(lun_postp) pleng
          read(lun_postp)  ( bridge(ileng), ileng = 1,pdime*pleng )

       end if

    else if( IMASTER ) then

       !-----------------------------------------------------------------
       !
       ! MASTER
       !
       !-----------------------------------------------------------------

       do kfl_desti_par = 1,npart

          if( kfl_reawr /= 1 ) then
             !
             ! - Master writes postpro without filter
             ! - Master writes restart
             !
             if(      wopos(3) == 'NPOIN' ) then
                if( associated(meshe) .and. size(meshe,kind=ip) >= imesh ) then
                   pleng = meshe(imesh) % npoin_par(kfl_desti_par)
                else
                   pleng = npoin_par(kfl_desti_par)
                end if

             else if( wopos(3) == 'NELEM' ) then
                pleng = meshe(imesh) % nelem_par(kfl_desti_par)
             else if( wopos(3) == 'NBOUN' ) then
                pleng = meshe(imesh) % nboun_par(kfl_desti_par)
             else if( wopos(3) == 'WHATE' ) then
                call PAR_SEND_RECEIVE(0_ip,1_ip,dummi,pleng_com,'IN MY CODE',kfl_desti_par)
                pleng = pleng_com(1)
             end if

             call postpr_memory(0_ip,max(1_ip,pdime*pleng),0_ip)

             call pararr('RCV',0_ip,pdime*pleng,gescp)

             write(lun_postp) pleng
             if( pleng > 0 ) write(lun_postp) ( gescp(ileng), ileng = 1,pdime*pleng )

             call postpr_memory(2_ip,max(1_ip,pdime*pleng),0_ip)

          else if( kfl_reawr == 1 ) then
             !
             ! - Master reads restart and sends to slave
             !
             if(      wopos(3) == 'NPOIN' ) then
                pleng = npoin_par(kfl_desti_par)
             else if( wopos(3) == 'NELEM' ) then
                pleng = nelem_par(kfl_desti_par)
             else if( wopos(3) == 'NBOUN' ) then
                pleng = nboun_par(kfl_desti_par)
             else if( wopos(3) == 'WHATE' ) then
                call PAR_SEND_RECEIVE(0_ip,1_ip,dummi,pleng_com,'IN MY CODE',kfl_desti_par)
                pleng = pleng_com(1)
             end if
             call postpr_memory(0_ip,max(1_ip,pdime*pleng),0_ip)
             read(lun_postp) pleng
             if( pleng > 0 ) read(lun_postp) ( gescp(ileng), ileng = 1,pdime*pleng )
             call pararr('SND',0_ip,pdime*pleng,gescp)
             call postpr_memory(2_ip,max(1_ip,pdime*pleng),0_ip)

          end if
       end do

    else if( ISLAVE ) then

       !-----------------------------------------------------------------
       !
       ! SLAVES
       !
       !-----------------------------------------------------------------

       kfl_desti_par = 0

       if( wopos(3) == 'WHATE' ) then
          pleng_com(1) = pleng
          call PAR_SEND_RECEIVE(1_ip,0_ip,pleng_com,dummi,'IN MY CODE',kfl_desti_par)
       end if

       if( kfl_reawr == 0 ) then
          !
          ! Slaves for postprocess
          !
          if(      kfl_permu == 1 ) then

             call postpr_memory(0_ip,pdime,pleng)
             call permut_rp2(pdime,pleng,lpmsh,bridge,gevep)
             call pararr('SND',0_ip,pdime*pleng,gevep)
             call postpr_memory(2_ip,pdime,pleng)

          else if( kfl_perme == 1 ) then

             call postpr_memory(0_ip,pdime,pleng)
             call permut_rp2(pdime,pleng,lemsh,bridge,gevep)
             call pararr('SND',0_ip,pdime*pleng,gevep)
             call postpr_memory(2_ip,pdime,pleng)

          else if( kfl_permb == 1 ) then

             call postpr_memory(0_ip,pdime,pleng)
             call permut_rp2(pdime,pleng,lbmsh,bridge,gevep)
             call pararr('SND',0_ip,pdime*pleng,gevep)
             call postpr_memory(2_ip,pdime,pleng)

          else
             
             call pararr('SND',0_ip,pdime*pleng,bridge)
             
          end if

       else if( kfl_reawr == 1 ) then
          !
          ! Slaves receive for restart
          !
          call pararr('RCV',0_ip,pdime*pleng,bridge)

       else if( kfl_reawr == 2 ) then
          !
          ! Slaves send for preliminary
          !
          call pararr('SND',0_ip,pdime*pleng,bridge)

       end if

    end if

    call postpr_close_file()

  end subroutine postpr_postprocess_rp

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Postprocess an integer array
  !> @details Postprocess an integer array
  !>
  !----------------------------------------------------------------------

  subroutine postpr_postprocess_ip(bridge,wopos,itste,ttime,pdime,TAG1,TAG2,mesh)

    character(5),                      intent(in)    :: wopos(*)
    integer(ip),                       intent(inout) :: bridge(*)
    integer(ip),                       intent(in)    :: itste
    real(rp),                          intent(in)    :: ttime
    integer(ip),                       intent(in)    :: pdime
    integer(ip),            optional,  intent(in)    :: TAG1
    integer(ip),            optional,  intent(in)    :: TAG2
    type(mesh_type_basic),  optional,  intent(in)    :: mesh
    integer(ip)                                      :: ileng,pleng
    integer(ip)                                      :: imesh,kinteger,tag1_opt,tag2_opt
    !
    ! Options
    !
    tag1_opt = 0
    tag2_opt = 0
    if( present(TAG1) ) tag1_opt = TAG1
    if( present(TAG2) ) tag2_opt = TAG2
    !
    ! Header
    !
    kinteger = kind(pdime)
    wwwww(1) = 'ALYA '
    wwwww(2) = 'V0003'
    wwwww(3) = wopos(1)
    wwwww(4) = wopos(2)
    wwwww(5) = wopos(3)
    iiiii(5) = int(tag1_opt,4)
    iiiii(6) = int(tag2_opt,4)
    wwwww(6) = 'INTEG'
    wwwww(7) = '4BYTE'
    iiiii(1) = int(pdime,4)
    iiiii(2) = 0_4
    iiiii(4) = int(itste,4)
    rrrrr(1) = ttime

    if( kinteger == 8_ip ) wwwww(7) = '8BYTE'
    wopos_pos(1) = wopos(1)
    !
    ! Mesh multiplication?
    !
    if( kfl_reawr == 0 ) then
       imesh = kfl_posdi
    else
       imesh = ndivi
    end if
    !
    ! Type
    !
    if( wopos(3) == 'NPOIN' ) then

       parii = NPOIN_TYPE
       pleng = meshe(imesh) % npoin
       if( IMASTER ) pleng = meshe(imesh) % npoin_total

    else if( wopos(3) == 'NELEM' ) then

       parii = NELEM_TYPE
       pleng = meshe(imesh) % nelem
       if( IMASTER ) pleng = meshe(imesh) % nelem_total

    else if( wopos(3) == 'NBOUN' ) then

       parii = NBOUN_TYPE
       pleng =  meshe(imesh) % nboun
       if( IMASTER ) pleng = meshe(imesh) % nboun_total

    else

       call runend('MOD_POSTPR - POSTPR_POSTPROCESS_IP: UNDEFINED POSTPROCESS')

    end if
    iiiii(2) = int(pleng,4)
    !
    ! Output solution
    !
    call postpr_open_file()
    if( kfl_reawr < 0 ) then
       kfl_reawr = abs(kfl_reawr)
       return
    end if
    !
    ! Read or write header
    !
    call postpr_header()

    if( ISEQUEN ) then

       if( kfl_reawr /= 1 ) then
          write(lun_postp) pleng
          write(lun_postp) ( bridge(ileng), ileng = 1,pdime*pleng )
       else if( kfl_reawr == 1 ) then
          read(lun_postp) pleng
          read(lun_postp)  ( bridge(ileng), ileng = 1,pdime*pleng )
       end if

    else if( IMASTER ) then

       do kfl_desti_par = 1,npart
          !
          ! Master 
          !
          if(      wopos(3) == 'NPOIN' ) then
             pleng = meshe(imesh) % npoin_par(kfl_desti_par)
          else if( wopos(3) == 'NELEM' ) then
             pleng = meshe(imesh) % nelem_par(kfl_desti_par)
          else if( wopos(3) == 'NBOUN' ) then
             pleng = meshe(imesh) % nboun_par(kfl_desti_par)
          end if
          call postpr_memory(1_ip,pdime*pleng,0_ip)
          if( pleng > 0 ) call parari('RCV',0_ip,pdime*pleng,giscp)
          write(lun_postp) pleng
          if( pleng > 0 ) write(lun_postp) ( giscp(ileng), ileng = 1,pdime*pleng )
          call postpr_memory(3_ip,pdime*pleng,0_ip)
       end do
       
    else if( ISLAVE ) then

       kfl_desti_par = 0
       !
       ! Slaves
       !
       if( pleng > 0 ) call parari('SND',0_ip,pdime*pleng,bridge)

    end if

    call postpr_close_file()

  end subroutine postpr_postprocess_ip

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-19
  !> @brief   Allocate and deallocate memory for generic scalar and vector arrays
  !> @details Allocate and deallocate memory for generic scalar and vector arrays:
  !>          ITASK=0 ... Allocate memory
  !>          ITASK=2 ... Deallocate memory
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_memory(itask,ndim1,ndim2)

    use def_parame
    use def_master
    use mod_memory, only : memory_alloca
    use mod_memory, only : memory_deallo
    implicit none
    integer(ip), intent(in) :: itask,ndim1,ndim2
!    integer(4)              :: istat

    select case(itask)

    case(0_ip)
       !
       ! Allocate memory for real
       !
       if(ndim1>0.and.ndim2==0) then
          call memory_alloca(memke,'GESCP','postpr_memory',gescp,ndim1)
       else if(ndim1>0.and.ndim2>0) then
          call memory_alloca(memke,'GEVEP','postpr_memory',gevep,ndim1,ndim2)
       else if(ndim1<0.and.ndim2>0) then
          call runend('POSTPR_MEMORY: NOT CODED')
          !allocate(getep(-ndim1,-ndim1,ndim2),stat=istat)
          !call memchk(zero,istat,memke,'GETEP','postpr_memory',getep)
       end if

    case(1_ip)
       !
       ! Allocate memory for integer
       !
       if(ndim1/=0.and.ndim2==0) then
          call memory_alloca(memke,'GISCP','postpr_memory',giscp,ndim1)
       else if(ndim1/=0.and.ndim2/=0) then
          call memory_alloca(memke,'GIVEP','postpr_memory',givep,ndim1,ndim2)
       end if

    case(2_ip)
       !
       ! Deallocate memory for real
       !
       if(ndim1>0.and.ndim2==0) then
          call memory_deallo(memke,'GESCP','postpr_memory',gescp)
       else if(ndim1>0.and.ndim2>0) then
          call memory_deallo(memke,'GEVEP','postpr_memory',gevep)
       else if(ndim1<0.and.ndim2>0) then
          call memory_deallo(memke,'GETEP','postpr_memory',getep)
       end if

    case(3_ip)
       !
       ! Deallocate memory for integer
       !
       if(ndim1/=0.and.ndim2==0) then
          call memory_deallo(memke,'GISCP','postpr_memory',giscp)
       else if(ndim1/=0.and.ndim2/=0) then
          call memory_deallo(memke,'GIVEP','postpr_memory',givep)
       end if

    case(4_ip)
       !
       ! Allocate memory for r3p
       !
       call runend('POSTPR_MEMORY: NOT PRGRAMMED')

    case(5_ip)
       !
       ! Deallocate memory for r3p
       !
       call runend('POSTPR_MEMORY: NOT PRGRAMMED')

    end select

  end subroutine postpr_memory

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-19
  !> @brief   Open file
  !> @details Open postprocess/restart file
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_open_file(TAG1,TAG2)

    use mod_iofile
    use mod_opfpos

    integer(ip),  intent(in), optional  :: TAG1
    integer(ip),  intent(in), optional  :: TAG2
    logical(lg)                         :: mesh
    character(200)                      :: messa
    
    call opfpos_name(wopos_pos(1), ".alyabin", mesh, TAG1, TAG2)
    
    if( mesh ) then
       !
       ! This is a mesh file
       !
       if(INOTSLAVE) then
          call iofile(0_ip,lun_postp,fil_postp,'ALYA POSTPROCESS FILE','unknown','unformatted')
       end if
       
    else if( kfl_reawr == 0 ) then
       !
       ! Open postprocess file name
       !
       if( INOTSLAVE ) then
          call iofile(0_ip,lun_postp,fil_postp,'ALYA POSTPROCESS FILE','unknown','unformatted')
          write(lun_pos01,'(a)') trim(fil_postp)
          call iofile_flush_unit(lun_pos01)
       end if

    else if( kfl_reawr == 1 ) then
       !
       ! Open restart file name for reading
       !
       if( INOTSLAVE )  then
          call iofile(4_ip,lun_postp,fil_postp,'ALYA RESTART FILE','unknown','unformatted')
       end if
       call PAR_BROADCAST(kfl_reawr)
       if( kfl_reawr < 0 ) then
          messa = 'CANNOT OPEN RESTART FILE: '//trim(fil_postp)
          call outfor(2_ip,0_ip,trim(messa))
          call livinf(0_ip,trim(messa),0_ip)
       end if
       if( kfl_reawr == 1 .and. INOTSLAVE ) then
          call iofile(0_ip,lun_postp,fil_postp,'ALYA RESTART FILE','unknown','unformatted')
       end if

    else if( kfl_reawr == 2 .and. INOTSLAVE ) then
       !
       ! Open restart file name for writing
       !
       call iofile(0_ip,lun_postp,fil_postp,'ALYA RESTART FILE','unknown','unformatted')
    end if

  end subroutine postpr_open_file

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-02-19
  !> @brief   Close file
  !> @details Close postprocess/restart file
  !>
  !-----------------------------------------------------------------------

  subroutine postpr_close_file()

    use mod_iofile

    if( INOTSLAVE ) then
       if( kfl_reawr == 0 ) then
          call iofile(2_ip,lun_postp,' ','ALYA POSTPROCESS FILE')
       else
          call iofile(2_ip,lun_postp,' ','ALYA RESTART FILE')
       end if
    end if

  end subroutine postpr_close_file

  subroutine postpr_read_restart_RP_1(bridge,ivari,itste,ttime,worig,pleng_opt,TAG1,TAG2,mesh)

    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    type(mesh_type_basic),  optional,         intent(in)    :: mesh

    kfl_reawr = 1
    call postpr_options(ivari)
    call postpr_real_scalar(bridge,postp(1) % wopos(:,ivari),itste,ttime,worig,pleng_opt,TAG1,TAG2,MESH=mesh)
    kfl_reawr = 0

  end subroutine postpr_read_restart_RP_1

  subroutine postpr_read_restart_IP_1(bridge,ivari,itste,ttime,worig,TAG1,TAG2)

    implicit none
    integer(ip),  intent(in)             :: ivari
    integer(ip),  intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2

    kfl_reawr = 1
    call postpr_options(ivari)
    call postpr_int_scalar(bridge,postp(1) % wopos(:,ivari),itste,ttime,worig,TAG1,TAG2)
    kfl_reawr = 0

  end subroutine postpr_read_restart_IP_1

  subroutine postpr_read_restart_IP_2(bridge,ivari,itste,ttime,kdime,worig,TAG1,TAG2)

    implicit none
    integer(ip),  intent(in)             :: ivari
    integer(ip),  intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2

    kfl_reawr = 1
    call postpr_options(ivari)
    call postpr_int_vector(bridge,postp(1) % wopos(:,ivari),itste,ttime,kdime,worig,TAG1,TAG2) 
    kfl_reawr = 0

  end subroutine postpr_read_restart_IP_2

  subroutine postpr_write_restart_RP_1(bridge,ivari,itste,ttime,worig,pleng_opt,TAG1,TAG2)

    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2

    kfl_reawr = 2
    call postpr_options(ivari)
    call postpr_real_scalar(bridge,postp(1) % wopos(:,ivari),itste,ttime,worig,pleng_opt,TAG1,TAG2)
    kfl_reawr = 0

  end subroutine postpr_write_restart_RP_1

  subroutine postpr_write_restart_IP_1(bridge,ivari,itste,ttime,worig,TAG1,TAG2)

    implicit none
    integer(ip),  intent(in)             :: ivari
    integer(ip),  intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2

    kfl_reawr = 2
    call postpr_options(ivari)
    call postpr_int_scalar(bridge,postp(1) % wopos(:,ivari),itste,ttime,worig,TAG1,TAG2)
    kfl_reawr = 0

  end subroutine postpr_write_restart_IP_1

  subroutine postpr_write_restart_IP_2(bridge,ivari,itste,ttime,kdime,worig,TAG1,TAG2)

    implicit none
    integer(ip),  intent(in)             :: ivari
    integer(ip),  intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime 
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2

    kfl_reawr = 2
    call postpr_options(ivari)
    call postpr_int_vector(bridge,postp(1) % wopos(:,ivari),itste,ttime,kdime,worig,TAG1,TAG2) 
    kfl_reawr = 0

  end subroutine postpr_write_restart_IP_2

  subroutine postpr_options(ivari,posit,posic,mesh)

    integer(ip),                     intent(in) :: ivari
    integer(ip),           optional, intent(in) :: posit
    integer(ip),           optional, intent(in) :: posic
    type(mesh_type_basic), optional, intent(in) :: mesh

    modul_loc      = modul
    time_posit_loc = momod(modul_loc) % postp(1) % time_posit(ivari)
    comp_posit_loc = momod(modul_loc) % postp(1) % comp_posit(ivari)
    enti_posit_loc = momod(modul_loc) % postp(1) % enti_posit(ivari)
    wopos_loc      = momod(modul_loc) % postp(1) % wopos(:,ivari)

    if(      momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NPOIN' ) then
       nenti = npoin
    else if( momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NELEM' ) then
       nenti = nelem
    else if( momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NBOUN' ) then
       nenti = nboun
    else
       call runend('MOD_POSTPR: UNKNOWN ENTITY TYPE '//trim(momod(modul_loc) % postp(1) % wopos(3,ivari)))
    end if

    if( present(mesh) ) then
       id_loc = mesh % id
    else
       id_loc = 1
    end if

  end subroutine postpr_options

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-14
  !> @brief   Read/write rp 
  !> @details Read write real(rp) :: xx(:,:)
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_read_restart_RP_2(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)

    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip),  intent(in), optional   :: posit

    kfl_reawr = 1
    call postpr_options(ivari,posit)
    call postpr_read_write_restart_RP_2(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)
    kfl_reawr = 0

  end subroutine postpr_read_restart_RP_2

  subroutine postpr_write_restart_RP_2(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)

    implicit none
    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip),  intent(in), optional   :: posit

    kfl_reawr = 2
    call postpr_options(ivari,posit)
    call postpr_read_write_restart_RP_2(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)
    kfl_reawr = 0

  end subroutine postpr_write_restart_RP_2

  subroutine postpr_read_write_restart_RP_2(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)

    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip),  intent(in), optional   :: posit
    integer(ip)                          :: itags,ienti,pdime
    integer(ip)                          :: itags_last,idime,ncomp,itags_first
    character(len=fsize_pos)             :: wfile
    logical(lg)                          :: read_write

    itags_last = 0
    
    if(      wopos_loc(2) == 'MATRI' ) then
       !
       ! Matrix
       !
       if( enti_posit_loc == 2 ) then
          pdime = memory_size(bridge,1_ip)
          call PAR_MAX(pdime)
          call postpr_real_vector(bridge,wopos_loc,itste,ttime,pdime,worig,pleng_opt)
       else
          call runend('postpr_read_write_restart_RP_2: DO NOT KNOW WHAT TO DO FOR '//wopos_loc(1))
       end if

    else if( wopos_loc(2) == 'VECTO' ) then
       !
       ! vector. Example: COORD(NDIME,NPOIN)
       !
       if( enti_posit_loc == 2 ) then       
          call postpr_real_vector(bridge,wopos_loc,itste,ttime,kdime,worig,pleng_opt)
       else
          call runend('postpr_read_write_restart_RP_2: DO NOT KNOW WHAT TO DO FOR '//wopos_loc(1))
       end if

    else if( wopos_loc(2) == 'SCALA' ) then
       !
       ! Scalar. Example: DEPOE(NTYLA,NELEM)
       !       
       if( enti_posit_loc == 2 .and. time_posit_loc == 0_ip ) then

          call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar,nenti)

          ncomp = memory_size(bridge,1_ip)
          call PAR_MAX(ncomp)

          do idime = 1,ncomp
             if( kfl_reawr == 2 ) then
                do ienti = 1,nenti
                   bridge_rp_scalar(ienti) = bridge(idime,ienti)
                end do
             end if
             if( present(TAG1) ) then
                call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=TAG1)
             else
                call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=idime)
             end if
             if( kfl_reawr == 1 ) then
                do ienti = 1,nenti
                   bridge(idime,ienti) = bridge_rp_scalar(ienti) 
                end do
             end if
          end do

          call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar)       

       else if( enti_posit_loc == 1 .and. time_posit_loc == 2_ip ) then
          !
          ! Scalar. Example: PRESS(NPOIN,NCOMP)
          !       
          call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar,nenti) 

          pdime = postp(1) % time_num(ivari)
          itags_first = 0
          
          do itags = 1,pdime

             read_write = .true.

             if( momod(modul_loc) % postp(1) % rst_time(min(itags,nvati),ivari) ) then
                if( itags_first == 0 ) itags_first = itags
                read_write = .true.
                !
                ! Copy to postprocess array
                !
                if( kfl_reawr == 2 ) then
                   do ienti = 1,nenti
                      bridge_rp_scalar(ienti) = bridge(ienti,itags)
                   end do
                end if
                !
                ! Check if file exists
                !
                if( kfl_reawr == 1 ) then
                   if( present(TAG1) ) then
                      call postpr_file_name(wfile,wopos_loc,TAG1=TAG1,TAG2=itags)
                   else
                      call postpr_file_name(wfile,wopos_loc,TAG1=1_ip,TAG2=itags)
                   end if
                   if( iofile_file_exists(wfile) ) then
                      itags_last = itags
                   else
                      call messages_live('RESTART FILE DOES NOT EXIST: '//trim(wfile)//'... COPYING LAST AVAILABLE TIME STEP','WARNING')
                      read_write = .false.
                   end if
                end if
                !
                ! Read or write
                !
                if( read_write ) then
                   if( present(TAG1) ) then
                      call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=TAG1,TAG2=itags)
                   else
                      call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=1_ip,TAG2=itags)
                   end if
                end if
                !
                ! Copy back if file has been read
                !
                if( kfl_reawr == 1 ) then
                   if( read_write ) then
                      do ienti = 1,nenti
                         bridge(ienti,itags) = bridge_rp_scalar(ienti) 
                      end do
                   else
                      if( itags_last == 0 ) call runend('MOD_POSTPR: WRONG ITAG_LAST')
                      do ienti = 1,nenti
                         bridge(ienti,itags) = bridge(ienti,itags_last) 
                      end do
                   end if
                end if

             else if( kfl_reawr == 1 ) then
                
                if( itags_first == 0 ) call runend('MOD_POSTPR: NO INFO TO COPY RESTART')
                call messages_live('COPY '//wopos_loc(1)//' (RST - PARALLEL) '//trim(postpr_tags(1_ip,itags))//' USING TIME STEP '//trim(intost(itags_first)))
                do ienti = 1,nenti
                   bridge(ienti,itags) = bridge(ienti,itags_first)
                end do

             end if

          end do
          call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar)       
       else
          call runend('MOD_POSTPR: WE ARE IN TROUBLE')
       end if
    else
       call runend('MOD_POSTPR: WE ARE IN TROUBLE')
    end if

  end subroutine postpr_read_write_restart_RP_2

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-14
  !> @brief   Read/write rp 
  !> @details Read write real(rp) :: xx(:,:,:)
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_read_restart_RP_3(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)

    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:,:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip),  intent(in), optional   :: posit

    kfl_reawr = 1
    call postpr_options(ivari,posit)
    call postpr_read_write_restart_RP_3(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)
    kfl_reawr = 0

  end subroutine postpr_read_restart_RP_3

  subroutine postpr_write_restart_RP_3(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)

    implicit none
    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:,:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip),  intent(in), optional   :: posit

    kfl_reawr = 2
    call postpr_options(ivari,posit)
    call postpr_read_write_restart_RP_3(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)
    kfl_reawr = 0

  end subroutine postpr_write_restart_RP_3

  subroutine postpr_read_write_restart_RP_3(bridge,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)

    integer(ip),  intent(in)             :: ivari
    real(rp),     intent(inout), pointer :: bridge(:,:,:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    integer(ip),  optional               :: kdime
    character(5), optional               :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip),  intent(in), optional   :: posit
    integer(ip)                          :: itags,idime,jdime,ienti,dim_time,ii
    integer(ip)                          :: ndim1,ndim3,ndim13
    integer(ip)                          :: itags_last,ncomp
    integer(ip)                          :: itags_first
    character(len=fsize_pos)             :: wfile
    logical(lg)                          :: read_write

    dim_time    = postp(1) % time_num(ivari)
    itags_first = 0          
    itags_last  = 0
    
    if( wopos_loc(2) == 'MATRI' ) then
       !
       ! Matrix
       !
       if( enti_posit_loc == 3 .and. time_posit_loc == 0_ip ) then
          !
          ! Does it really exist?
          !
          call postpr_real_matrix(bridge,wopos_loc,itste,ttime,kdime,worig,pleng_opt)

       else if( enti_posit_loc == 2 .and. time_posit_loc == 0_ip ) then
          ndim1 = memory_size(bridge,1_ip)
          ndim3 = memory_size(bridge,3_ip)
          ndim13 = ndim1*ndim3
          call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector,ndim13,nenti)
          if( kfl_reawr == 2 ) then
             do ienti = 1,nenti
                ii = 0
                do idime = 1,ndim3
                   do jdime = 1,ndim1
                      ii = ii + 1
                      bridge_rp_vector(ii,ienti) = bridge(jdime,ienti,idime)
                   end do
                end do
             end do
          end if
          call postpr_real_vector(bridge_rp_vector,wopos_loc,itste,ttime,ndim13,worig,pleng_opt)

          if( kfl_reawr == 1 ) then            
             do ienti = 1,nenti
                ii = 0
                do idime = 1,ndim3
                   do jdime = 1,ndim1
                      ii = ii + 1
                      bridge(jdime,ienti,idime) = bridge_rp_vector(ii,ienti)
                   end do
                end do
             end do
          end if
          call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector)

       else if( enti_posit_loc == 1 .and. time_posit_loc == 3_ip ) then
          !
          ! Matrix: XXXX(NPOIN,MGAUS*NDIME,NCOMP)
          !
          if( present(kdime) ) then
             call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector,kdime,nenti)

             do itags = 1,dim_time
                if( kfl_reawr == 2 ) then
                   do ienti = 1,nenti
                      do idime = 1,kdime
                         bridge_rp_vector(idime,ienti) = bridge(ienti,idime,itags)
                      end do
                   end do
                end if
                call postpr_real_vector(bridge_rp_vector,wopos_loc,itste,ttime,kdime,worig,pleng_opt,TAG1=1_ip,TAG2=itags)
                if( kfl_reawr == 1 ) then
                   do ienti = 1,nenti
                      do idime = 1,kdime
                         bridge(ienti,idime,itags) = bridge_rp_vector(idime,ienti) 
                      end do
                   end do
                end if
             end do
             call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector)
          else
             call runend('POSTPR: KDIME NOT PRESENT')
          end if

       else if( enti_posit_loc == 2 .and. time_posit_loc == 3_ip ) then
          !
          ! Matrix: XXXX(MGAUS*NDIME,NPOIN,NCOMP)
          !
          if( present(kdime) ) then

             call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector,kdime,nenti)

             do itags = 1,dim_time
                if( kfl_reawr == 2 ) then
                   do ienti = 1,nenti
                      do idime = 1,kdime
                         bridge_rp_vector(idime,ienti) = bridge(idime,ienti,itags)
                      end do
                   end do
                end if
                call postpr_real_vector(bridge_rp_vector,wopos_loc,itste,ttime,kdime,worig,pleng_opt,TAG1=1_ip,TAG2=itags)
                if( kfl_reawr == 1 ) then
                   do ienti = 1,nenti
                      do idime = 1,kdime
                         bridge(idime,ienti,itags) = bridge_rp_vector(idime,ienti) 
                      end do
                   end do
                end if
             end do
             call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector)
          else
             call runend('POSTPR: KDIME NOT PRESENT')
          end if

       else
          call runend('MOD_POSTPR: WE ARE IN TROUBLE WITH A MATRIX')          
       end if

    else if( wopos_loc(2) == 'VECTO' ) then 
       !
       ! Vector. Example: VELOC(NDIME,NPOIN,NCOMP)
       !
       if( enti_posit_loc == 2 .and. time_posit_loc == 3_ip ) then

          call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector,ndime,nenti)

          do itags = 1,dim_time

             if( momod(modul_loc) % postp(1) % rst_time(min(itags,nvati),ivari) ) then
                if( itags_first == 0 ) itags_first = itags
                read_write = .true.
                !
                ! Copy to postprocess array
                !           
                if( kfl_reawr == 2 ) then
                   do ienti = 1,nenti
                      do idime = 1,memory_size(bridge,1_ip)
                         bridge_rp_vector(idime,ienti) = bridge(idime,ienti,itags)
                      end do
                   end do
                end if
                !
                ! Check if file exists
                !
               
                if( kfl_reawr == 1 ) then
                   if( present(TAG1) ) then
                      call postpr_file_name(wfile,wopos_loc,TAG1=TAG1,TAG2=itags)
                   else
                      call postpr_file_name(wfile,wopos_loc,TAG1=1_ip,TAG2=itags)
                   end if
                   if( iofile_file_exists(wfile) ) then
                      itags_last = itags
                   else
                      call messages_live('RESTART FILE DOES NOT EXIST: '//trim(wfile)//'... COPYING LAST AVAILABLE TIME STEP','WARNING')
                      read_write = .false.
                   end if
                end if                
                !
                ! Read or write
                !             
                if( read_write ) then
                   if( present(TAG1) ) then
                      call postpr_real_vector(bridge_rp_vector,wopos_loc,itste,ttime,kdime,worig,pleng_opt,TAG1=TAG1,TAG2=itags)
                   else
                      call postpr_real_vector(bridge_rp_vector,wopos_loc,itste,ttime,kdime,worig,pleng_opt,TAG1=1_ip,TAG2=itags)
                   end if
                end if
                !
                ! Copy back if file has been read
                !            
                if( kfl_reawr == 1 ) then
                   if( read_write ) then
                      do ienti = 1,nenti
                         do idime = 1,memory_size(bridge,1_ip)
                            bridge(idime,ienti,itags) = bridge_rp_vector(idime,ienti) 
                         end do
                      end do
                   else
                      if( itags_last == 0 ) call runend('MOD_POSTPR: DID NOT FIND ANYTHING TO COPY')
                      do ienti = 1,nenti
                         do idime = 1,memory_size(bridge,1_ip)
                            bridge(idime,ienti,itags) = bridge(idime,ienti,itags_last)
                         end do
                      end do
                   end if
                end if

             else if( kfl_reawr == 1 ) then
                
                if( itags_first == 0 ) call runend('MOD_POSTPR: NO INFO TO COPY RESTART')
                call messages_live('COPY '//wopos_loc(1)//' (RST - PARALLEL) '//trim(postpr_tags(1_ip,itags))//' USING TIME STEP '//trim(intost(itags_first)))
                do ienti = 1,nenti
                   do idime = 1,memory_size(bridge,1_ip)
                      bridge(idime,ienti,itags) = bridge(idime,ienti,itags_first)
                   end do
                end do
                
             end if

          end do
         
          call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_vector)       

       else

          call runend('MOD_POSTPR: WE ARE IN TROUBLE WITH A MATRIX '//wopos_loc(1))

       end if

    else if( wopos_loc(2) == 'SCALA' ) then

       if( enti_posit_loc == 1 .and. time_posit_loc == 3_ip ) then
          !
          ! Scalar. Example: CONCE(NPOIN,NCLAS,NCOMP) 
          !
          call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar,nenti)
          ncomp = memory_size(bridge,2_ip)
          call PAR_MAX(ncomp)

          do itags = 1,dim_time
             
             if( momod(modul_loc) % postp(1) % rst_time(min(itags,nvati),ivari) ) then
                if( itags_first == 0 ) itags_first = itags
                do idime = 1,ncomp
                   if( kfl_reawr == 2 ) then
                      do ienti = 1,nenti
                         bridge_rp_scalar(ienti) = bridge(ienti,idime,itags)
                      end do
                   end if
                   if( present(TAG1) ) then
                      call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=TAG1,TAG2=itags)
                   else
                      call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=idime,TAG2=itags)
                   end if
                   if( kfl_reawr == 1 ) then
                      do ienti = 1,nenti
                         bridge(ienti,idime,itags) = bridge_rp_scalar(ienti) 
                      end do
                   end if
                end do

             else if( kfl_reawr == 1 ) then
                
                if( itags_first == 0 ) call runend('MOD_POSTPR: NO INFO TO COPY RESTART')
                do idime = 1,ncomp
                   call messages_live('COPY '//wopos_loc(1)//' (RST - PARALLEL) '//trim(postpr_tags(idime,itags))//' USING TIME STEP '//trim(intost(itags_first)))
                   do ienti = 1,nenti
                      bridge(ienti,idime,itags) = bridge(ienti,idime,itags_first)
                   end do
                end do
                
             end if
             
          end do

          call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar)       

       else if( enti_posit_loc == 2 .and. time_posit_loc == 3_ip ) then
          !
          ! Scalar. Example: VAUXI_EXM(NAUXI_EXM,NPOIN,NCOMP) 
          !
          call memory_alloca(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar,nenti)
          ncomp = memory_size(bridge,1_ip)
          call PAR_MAX(ncomp)

          do itags = 1,dim_time
             
             if( momod(modul_loc) % postp(1) % rst_time(min(itags,nvati),ivari) ) then
                if( itags_first == 0 ) itags_first = itags

                do idime = 1,ncomp
                   if( kfl_reawr == 2 ) then
                      do ienti = 1,nenti
                         bridge_rp_scalar(ienti) = bridge(idime,ienti,itags)
                      end do
                   end if
                   if( present(TAG1) ) then
                      call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=TAG1,TAG2=itags)
                   else
                      call postpr_real_scalar(bridge_rp_scalar,wopos_loc,itste,ttime,worig,pleng_opt,TAG1=idime,TAG2=itags)
                   end if
                   if( kfl_reawr == 1 ) then
                      do ienti = 1,nenti
                         bridge(idime,ienti,itags) = bridge_rp_scalar(ienti) 
                      end do
                   end if
                end do
                
             else if( kfl_reawr == 1 ) then
                
                if( itags_first == 0 ) call runend('MOD_POSTPR: NO INFO TO COPY RESTART')
                do idime = 1,ncomp
                   call messages_live('COPY '//wopos_loc(1)//' (RST - PARALLEL) '//trim(postpr_tags(idime,itags))//' USING TIME STEP '//trim(intost(itags_first)))
                   do ienti = 1,nenti
                      bridge(idime,ienti,itags) = bridge(idime,ienti,itags_first)
                   end do
                end do
                
             end if
             
          end do

          call memory_deallo(memke,wopos_loc(1),'mod_postpr',bridge_rp_scalar)       

       else
          call runend('MOD_POSTPR: WE ARE IN TROUBLE')
       end if

    else

       call runend('MOD_POSTPR: WE ARE IN TROUBLE WITH A MATRIX '//wopos_loc(1))

    end if

  end subroutine postpr_read_write_restart_RP_3

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-14
  !> @brief   Read/write r3p 
  !> @details Read write type(r3p) :: xx(:)
  !> 
  !-----------------------------------------------------------------------

  subroutine postpr_read_restart_R3P_1(bridge,ivari,itste,ttime,worig,pleng_opt,TAG1,TAG2)

    integer(ip),  intent(in)             :: ivari
    type(r3p),    intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5),             optional   :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    
    kfl_reawr = 1
    call postpr_options(ivari)
    call postpr_read_write_restart_R3P_1(bridge,ivari,itste,ttime,worig,pleng_opt,TAG1,TAG2)
    kfl_reawr = 0

  end subroutine postpr_read_restart_R3P_1

  subroutine postpr_write_restart_R3P_1(bridge,ivari,itste,ttime,worig,pleng_opt,TAG1,TAG2)

    integer(ip),  intent(in)             :: ivari
    type(r3p),    intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5),             optional   :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2

    kfl_reawr = 2
    call postpr_options(ivari)
    call postpr_read_write_restart_R3P_1(bridge,ivari,itste,ttime,worig,pleng_opt,TAG1,TAG2)
    kfl_reawr = 0

  end subroutine postpr_write_restart_R3P_1

  subroutine postpr_read_write_restart_R3P_1(bridge,ivari,itste,ttime,worig,pleng_opt,TAG1,TAG2,posit)

    integer(ip),  intent(in)             :: ivari
    type(r3p),    intent(inout), pointer :: bridge(:)
    integer(ip) , intent(in)             :: itste
    real(rp)    , intent(in)             :: ttime
    character(5),             optional   :: worig
    integer(ip),  intent(in), optional   :: pleng_opt
    integer(ip),  intent(in), optional   :: TAG1
    integer(ip),  intent(in), optional   :: TAG2
    integer(ip),  intent(in), optional   :: posit
    integer(ip)                          :: dim_time,kdime,idime,ii,jdime
    integer(ip)                          :: idim1,idim2,idim3,ndim1,ndim2
    real(rp),                 pointer    :: xx(:,:,:)
    
    nullify(xx)

    if( time_posit_loc == 3 ) then
       dim_time = 1
       kdime    = 0
       if( associated(bridge) ) then
          if( associated(bridge(1) % a) ) then
             dim_time = size(bridge(1) % a,dim=3)
             do ii = 1,nenti
                idime = size(bridge(ii) % a,DIM=1,KIND=ip)*size(bridge(ii) % a,DIM=2,KIND=ip)
                kdime = max(kdime,idime)
             end do
          end if
       end if       
       call PAR_MAX(dim_time)
       call PAR_MAX(kdime)
       wopos_loc(2)               = 'MATRI'
       enti_posit_loc             = 2
       postp(1) % time_num(ivari) = dim_time
       
       call memory_alloca(memke,wopos_loc(1),'mod_postpr',xx,kdime,nenti,dim_time)
       !
       ! Write
       !       
       if( kfl_reawr == 2 ) then
          do idim3 = 1,dim_time
             do ii = 1,nenti
                ndim2 = size(bridge(ii) % a,2)
                do idim2 = 1,ndim2
                   ndim1 = size(bridge(ii) % a,1)
                   do idim1 = 1,ndim1
                      jdime = (idim2-1)*ndim1+idim1
                      xx(jdime,ii,idim3) = bridge(ii) % a(idim1,idim2,idim3) 
                   end do
                end do
             end do
          end do
       end if

       call postpr_read_write_restart_RP_3(xx,ivari,itste,ttime,kdime,worig,pleng_opt,TAG1,TAG2,posit)
       !
       ! Read
       !
       if( kfl_reawr == 1 ) then
          do idim3 = 1,dim_time
             do ii = 1,nenti
                ndim2 = size(bridge(ii) % a,2)
                do idim2 = 1,ndim2
                   ndim1 = size(bridge(ii) % a,1)
                   do idim1 = 1,ndim1
                      jdime = (idim2-1)*ndim1+idim1
                      bridge(ii) % a(idim1,idim2,idim3) = xx(jdime,ii,idim3) 
                   end do
                end do
             end do
          end do
       end if
       
       call memory_deallo(memke,wopos_loc(1),'mod_postpr',xx)
       
    else
       call runend('postpr_read_write_restart_R3P_1: DO NOT KNOW WHAT TO DO FOR '//wopos_loc(1))
    end if

  end subroutine postpr_read_write_restart_R3P_1

end module mod_postpr
!> @}

