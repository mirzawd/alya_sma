!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_postpr.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO post process (parallel)
!> @details This module writes/reads restart and post process files in parallel using the MPI-IO format
!-----------------------------------------------------------------------


module mod_mpio_par_postpr

  use def_kintyp,             only : ip,rp,lg
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use def_postpr
  use def_mpio
  use def_parall
  use mod_mpio_seq_postpr
  use mod_mpio_par_io
  use mod_mpio_seq_io
  use mod_mpio_par_hybrid_io
  use mod_postpr_tools,       only : postpr_tags
  use mod_parall,             only : PAR_CODE_SIZE, PAR_MY_CODE_RANK
  use mod_parall,             only : PAR_COMM_MPIO_WM_SIZE
  use mod_communications,     only : PAR_BROADCAST, PAR_BARRIER
  use mod_communications,     only : PAR_COMM_SPLIT, PAR_GATHERV, PAR_GATHER, PAR_SUM
  use mod_messages,           only : messages_live
  use mod_memory,             only : memory_alloca, memory_deallo
  use mod_redistribution,     only : redistribution_array
  use mod_io,                 only : io_read_merge
  use def_mpi
  use mod_mpio_config,        only : mpio_config
#ifdef ALYA_FTI
  use mod_alya2fti,           only : FTI_st, FTI_write_ckpt
  use FTI
#endif
  
  implicit none

  private


  character(150)                :: wherein_code="IN MY CODE"
  type(comm_data_par), pointer  :: comm_postpr
  type(comm_data_par), pointer  :: comm_npoin_postpr
  type(comm_data_par), pointer  :: comm_nelem_postpr
  type(comm_data_par), pointer  :: comm_nboun_postpr
  type(comm_data_par), target   :: comm_npoin_postpr_rst
  type(comm_data_par), target   :: comm_nelem_postpr_rst
  type(comm_data_par), target   :: comm_nboun_postpr_rst
  type(comm_data_par), target   :: comm_npoin_postpr_post
  type(comm_data_par), target   :: comm_nelem_postpr_post
  type(comm_data_par), target   :: comm_nboun_postpr_post
  integer(ip)                   :: chunk_npoin_post
  integer(ip)                   :: chunk_nelem_post
  integer(ip)                   :: chunk_nboun_post
  integer(ip)                   :: chunk_npoin_rst
  integer(ip)                   :: chunk_nelem_rst
  integer(ip)                   :: chunk_nboun_rst
  logical(lg)                   :: init_redistrib_rst=.false.
  logical(lg)                   :: init_redistrib_post=.false.
  character(50)                 :: wtag_mpio

  public :: par_posmpio_int_m, par_posmpio_int_v, par_posmpio_real_m, par_posmpio_real_v
  public :: init_redistrib_post
  public :: init_redistrib_rst
  public :: posmpio_redistribute
  public :: posmpio_destructor

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-04
  !> @brief   Tags fro MPIO
  !> @details Tags fro MPIO
  !> 
  !-----------------------------------------------------------------------

  subroutine par_posmpio_tags()

    wtag_mpio = postpr_tags(tag1_mpio,tag2_mpio)

  end subroutine par_posmpio_tags  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-04
  !> @brief   MPIO
  !> @details INTEGER(:)
  !> 
  !-----------------------------------------------------------------------

  subroutine par_posmpio_int_v()
    type(mpio_header)                   :: header
    logical                             :: too_small
    integer(ip)                         :: my_nsubd
    character(5)                        :: postrst
#ifdef ALYA_FTI
    integer                             :: error
    integer                             :: FTI_id
    integer(ip), target                 :: FTI_dummy
    integer(ip), pointer                :: FTI_ptr_dummy

    FTI_ptr_dummy => FTI_dummy
#endif
#ifndef MPI_OFF
    if (kfl_reawr==0) then
        postrst='POST'
    else 
        postrst='RST'
     end if

    call mpio_postpr_broadcast()
    call mpio_postpr_init()
    call par_posmpio_tags()
    too_small=is_too_small(dim, ip)
    my_nsubd=PAR_CODE_SIZE-1
    
    if (kfl_reawr/=1) then
       !
       ! Write
       !
       call mpio_postpr_permu_int_v()
       if (mpio_config%output%merge) then
          if (INOTMASTER) then
             call redistribute_lelbo()
             call posmpio_redistribute()
             call redistribution_array(isca,resultson(1:5),MEMOR=memor_dom,VARIABLE_NAME='isca',COMM=comm_postpr,PERMUTE_RECV=.true.)
          end if
          my_nsubd=1
       end if
#ifdef ALYA_FTI
       if ( FTI_write_ckpt == 1 ) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - FTI)'//' '//trim(wtag_mpio))
          call FTI_setIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, isca, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_AddVarICP(FTI_id, error)
       else &
#endif
       if (is_parallel() .and. .not.too_small) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - PARALLEL)'//' '//trim(wtag_mpio))
          call PAR_FILE_WRITE_ALL(isca, fil_postp, object, resultson, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       else
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call PARSEQ_FILE_WRITE_ALL(isca, fil_postp, object, resultson, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       end if
       call mpio_postpr_permu_end()
       
    else if (kfl_reawr==1) then
       !
       ! Read
       !
#ifdef ALYA_FTI
       if (FTI_st /= 0 ) then
          call messages_live('READ '//object(1:5)//' (RST - FTI)'//' '//trim(wtag_mpio))
          call FTI_getIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, gescai_mpio, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_RecoverVar(FTI_id, error)
       else &
#endif
       if (is_parallel() .and. .not.too_small) then
          
          call messages_live('READ '//object(1:5)//' (RST - PARALLEL)'//' '//trim(wtag_mpio))

          if( mpio_config%output%merge ) then
             call PAR_FILE_READ_HEADER(fil_postp, header, wherein_code)
             call mpio_postpr_header(header)
             call io_read_merge(gescai_mpio,fil_postp,resultson)
          else
             call PAR_FILE_READ_HEADER(fil_postp, header, wherein_code)
             call mpio_postpr_header(header)
             call PAR_FILE_READ_ALL(gescai_mpio, fil_postp, dim, wherein_code, header=header)
          end if
          
       else
          
          call messages_live('READ '//object(1:5)//' (RST - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call SEQ_FILE_READ_HEADER(fil_postp, header)
          call mpio_postpr_header(header)
          call PARSEQ_FILE_READ_ALL(gescai_mpio, fil_postp, dim, header=header)
          
       end if
       
    end if
#endif
  end subroutine par_posmpio_int_v

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-04
  !> @brief   MPIO
  !> @details INTEGER(:,:)
  !> 
  !-----------------------------------------------------------------------

  subroutine par_posmpio_int_m()
    type(mpio_header)                   :: header
    logical                             :: too_small
    integer(ip)                         :: my_nsubd
    character(5)                        :: postrst
#ifdef ALYA_FTI
    integer                             :: error
    integer                             :: FTI_id
    integer(ip), target                 :: FTI_dummy
    integer(ip), pointer                :: FTI_ptr_dummy

    FTI_ptr_dummy => FTI_dummy
#endif
#ifndef MPI_OFF
    if (kfl_reawr==0) then
        postrst='POST'
    else
        postrst='RST'
    end if
    call mpio_postpr_broadcast()
    call mpio_postpr_init()
    call par_posmpio_tags()
    too_small = is_too_small(dim*pdime_mpio, ip)
    my_nsubd  = PAR_CODE_SIZE-1
    
    if (kfl_reawr/=1) then
       call mpio_postpr_permu_int_m()
       if (mpio_config%output%merge) then
          if (INOTMASTER) then
             call redistribute_lnods()
             call redistribute_lnodb()
             call posmpio_redistribute()
             !call redistribution_array(ivec,resultson(1:5),MEMOR=memor_dom,VARIABLE_NAME=object(1:5),COMM=comm_postpr,PERMUTE_RECV=.true.)
             call redistribution_array(ivec,resultson(1:5),MEMOR=memor_dom,VARIABLE_NAME='ivec',COMM=comm_postpr,PERMUTE_RECV=.true.)
          end if
          my_nsubd=1
       end if
#ifdef ALYA_FTI
       if (FTI_write_ckpt == 1 ) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - FTI)'//' '//trim(wtag_mpio))
          call FTI_setIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, ivec, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_AddVarICP(FTI_id,error)
       else &
#endif
       if (is_parallel() .and. .not.too_small) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - PARALLEL)'//' '//trim(wtag_mpio))
          call PAR_FILE_WRITE_ALL(ivec, fil_postp, object, resultson, pdime_mpio, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       else
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call PARSEQ_FILE_WRITE_ALL(ivec, fil_postp, object, resultson, pdime_mpio, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       end if
       call mpio_postpr_permu_end()
    else if (kfl_reawr==1) then
#ifdef ALYA_FTI
       if ( FTI_st /= 0 ) then
          call messages_live('READ '//object(1:5)//' (RST - FTI)'//' '//trim(wtag_mpio))
          call FTI_getIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, geveci_mpio, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_RecoverVar(FTI_id, error)
       else &
#endif
       if (is_parallel() .and. .not.too_small) then
          call messages_live('READ '//object(1:5)//' (RST - PARALLEL)'//' '//trim(wtag_mpio))
          call PAR_FILE_READ_HEADER(fil_postp, header, wherein_code)
          call mpio_postpr_header(header)
          call PAR_FILE_READ_ALL(geveci_mpio, fil_postp, pdime_mpio, dim, wherein_code, header=header)
       else
          call messages_live('READ '//object(1:5)//' (RST - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call SEQ_FILE_READ_HEADER(fil_postp, header)
          call mpio_postpr_header(header)
          call PARSEQ_FILE_READ_ALL(geveci_mpio, fil_postp, pdime_mpio, dim, header=header)
       end if
    end if
#endif
  end subroutine par_posmpio_int_m
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-04
  !> @brief   MPIO
  !> @details REAL(:)
  !> 
  !-----------------------------------------------------------------------

  subroutine par_posmpio_real_v()
    type(mpio_header)                   :: header
    logical                             :: too_small
    integer(ip)                         :: my_nsubd
    character(5)                        :: postrst
#ifdef ALYA_FTI
    integer                             :: error
    integer                             :: FTI_id
    real(rp), target                    :: FTI_dummy
    real(rp), pointer                   :: FTI_ptr_dummy

    FTI_ptr_dummy => FTI_dummy
#endif

#ifndef MPI_OFF
    if (kfl_reawr==0) then
        postrst='POST'
    else
        postrst='RST'
    end if
    call mpio_postpr_broadcast()
    call mpio_postpr_init()
    call par_posmpio_tags()
    too_small = is_too_small(dim, rp)
    my_nsubd  = PAR_CODE_SIZE-1
    if (kfl_reawr/=1) then
       !
       ! Write post (kfl_reawr=0) or write restart (kfl_reawr=2)
       !
       call mpio_postpr_permu_real_v()
       if (mpio_config%output%merge) then
          if (INOTMASTER) then
             call posmpio_redistribute()
             call redistribution_array(rsca,resultson(1:5),MEMOR=memor_dom,VARIABLE_NAME='rsca',COMM=comm_postpr,PERMUTE_RECV=.true.)
          end if
          my_nsubd=1
       end if
#ifdef ALYA_FTI
       if ( FTI_write_ckpt == 1 ) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - FTI)'//' '//trim(wtag_mpio))
          call FTI_setIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, rsca, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_AddVarICP(FTI_id, error)
!          if ( dim /= 0 ) then
!            print *,rsca(1:5)
!          end if
       else &
#endif
       if (is_parallel() .and. .not.too_small) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - PARALLEL)'//' '//trim(wtag_mpio))
          call PAR_FILE_WRITE_ALL(rsca, fil_postp, object, resultson, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       else
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call PARSEQ_FILE_WRITE_ALL(rsca, fil_postp, object, resultson, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       end if
       call mpio_postpr_permu_end()
       
    else if (kfl_reawr==1) then
       !
       ! Read restart
       !       
#ifdef ALYA_FTI
       if ( FTI_st /= 0 ) then
          call messages_live('READ '//object(1:5)//' (RST - FTI)'//' '//trim(wtag_mpio))
          call FTI_getIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, gescar_mpio, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_RecoverVar(FTI_id, error)
!          if ( dim /= 0 ) then
!            print *,gescar_mpio(1:5)
!          end if
       else &
#endif
       if (is_parallel() .and. .not.too_small) then
          
          call messages_live('READ '//object(1:5)//' (RST - PARALLEL)'//' '//trim(wtag_mpio))

          if( mpio_config%output%merge ) then
             call PAR_FILE_READ_HEADER(fil_postp, header, wherein_code)
             call mpio_postpr_header(header)
             call io_read_merge(gescar_mpio,fil_postp,resultson)
          else
             call PAR_FILE_READ_HEADER(fil_postp, header, wherein_code)
             call mpio_postpr_header(header)
             call PAR_FILE_READ_ALL(gescar_mpio, fil_postp, dim, wherein_code, header=header)
          end if
       else
          
          call messages_live('READ '//object(1:5)//' (RST - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call SEQ_FILE_READ_HEADER(fil_postp, header)
          call mpio_postpr_header(header)
          call PARSEQ_FILE_READ_ALL(gescar_mpio, fil_postp, dim, header=header)
          
       end if
    end if
#endif
  end subroutine par_posmpio_real_v

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-04
  !> @brief   MPIO
  !> @details REAL(:,:)
  !> 
  !-----------------------------------------------------------------------

  subroutine par_posmpio_real_m()
    use mod_io
    type(mpio_header)                   :: header
    logical                             :: too_small
    integer(ip)                         :: my_nsubd
    character(5)                        :: postrst
#ifdef ALYA_FTI
    integer                             :: error
    integer                             :: FTI_id
    real(rp), target                    :: FTI_dummy
    real(rp), pointer                   :: FTI_ptr_dummy
    FTI_ptr_dummy => FTI_dummy
#endif
#ifndef MPI_OFF
    if( kfl_reawr == 0 ) then
        postrst='POST'
    else
        postrst='RST'
    end if
    call mpio_postpr_broadcast()
    call mpio_postpr_init()
    call par_posmpio_tags()
    too_small=is_too_small(dim*pdime_mpio, rp)

    my_nsubd = PAR_CODE_SIZE-1
    
    if( kfl_reawr /= 1 ) then
       !
       ! Write post (kfl_reawr=0) or write restart (kfl_reawr=2)
       !
       call mpio_postpr_permu_real_m()
       if( mpio_config%output%merge ) then
          if( INOTMASTER ) then
             call posmpio_redistribute()
             call redistribution_array(rvec,resultson(1:5),MEMOR=memor_dom,VARIABLE_NAME='rvec',COMM=comm_postpr,PERMUTE_RECV=.true.)
          end if
          my_nsubd=1
       end if
#ifdef ALYA_FTI
       if ( FTI_write_ckpt == 1 ) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - FTI)'//' '//trim(wtag_mpio))
          call FTI_setIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, rvec, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_AddVarICP(FTI_id, error)
       else &
#endif
       if (is_parallel() .and. .not.too_small) then
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - PARALLEL)'//' '//trim(wtag_mpio))
          call PAR_FILE_WRITE_ALL(rvec, fil_postp, object, resultson, pdime_mpio, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       else
          call messages_live('WRITE '//object(1:5)//' ('//trim(postrst)//' - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call PARSEQ_FILE_WRITE_ALL(rvec, fil_postp, object, resultson, pdime_mpio, dim, imesh, tag1_mpio, tag2_mpio, itste_mpio, ttime_mpio, nsubd=my_nsubd)
       end if
       call mpio_postpr_permu_end()
       
    else if( kfl_reawr == 1 ) then
       !
       ! Read restart
       !
#ifdef ALYA_FTI
       if ( FTI_st /= 0 ) then
          call messages_live('READ '//object(1:5)//' (RST - FTI)'//' '//trim(wtag_mpio))
          call FTI_getIDFromString(object(1:5), FTI_id)
          if ( dim /= 0 ) then
             call FTI_Protect(FTI_id, gevecr_mpio, error)
          else
             call FTI_Protect(FTI_id, FTI_ptr_dummy, error)
          end if
          call FTI_RecoverVar(FTI_id, error)
       else &
#endif
       if( is_parallel() .and. .not. too_small ) then
          
          call messages_live('READ '//object(1:5)//' (RST - PARALLEL)'//' '//trim(wtag_mpio))
          
          if( mpio_config%output%merge ) then
             call PAR_FILE_READ_HEADER(fil_postp, header, wherein_code)
             call mpio_postpr_header(header)
             call io_read_merge(gevecr_mpio,fil_postp,resultson)
          else
             call PAR_FILE_READ_HEADER(fil_postp, header, wherein_code)
             call mpio_postpr_header(header)
             call PAR_FILE_READ_ALL(gevecr_mpio, fil_postp, pdime_mpio, dim, wherein_code, header=header)
          end if
       else
          call messages_live('READ '//object(1:5)//' (RST - SEQUENTIAL (FORCED))'//' '//trim(wtag_mpio))
          call SEQ_FILE_READ_HEADER(fil_postp, header)
          call mpio_postpr_header(header)
          call PARSEQ_FILE_READ_ALL(gevecr_mpio, fil_postp, pdime_mpio, dim, header=header)
       end if
       
    end if
#endif
  end subroutine par_posmpio_real_m

  subroutine mpio_postpr_broadcast()
    real(rp), pointer    ::  t(:)
    nullify(t)
    call memory_alloca(mpio_memor,'t','mpio_postpr_broadcast',t,PAR_CODE_SIZE+1)
    call PAR_BROADCAST(5_ip, wopos_mpio(1), wherein_code)
    call PAR_BROADCAST(5_ip, wopos_mpio(2), wherein_code)
    call PAR_BROADCAST(5_ip, wopos_mpio(3), wherein_code)
    call PAR_GATHER(ttime_mpio, t, wherein_code)
    if (IMASTER) then
       ttime_mpio=t(2)
    end if
    call memory_deallo(mpio_memor,'t','mpio_postpr_broadcast',t)
  end subroutine mpio_postpr_broadcast

  function is_too_small(dims, type_size) result(too_small)
    integer(ip), intent(in)    ::  dims
    logical(lg)                ::  too_small
    integer(ip)                ::  dimstmp
    integer(ip)                ::  avgsize
    integer, intent(in)        ::  type_size
    if (mpio_config%output%par_min_block_size .le. 0_ip) then
       too_small=.false.
    else
       dimstmp=dims
       if (IMASTER) then
          dimstmp=0
       end if
       call PAR_SUM(dimstmp, wherein_code)
       avgsize=int((real(type_size,rp)*real(dimstmp,rp))/((1024.0_rp*1024.0_rp*(real(PAR_CODE_SIZE,rp)-1.0_rp))), ip)
       call PAR_BROADCAST(avgsize, wherein=wherein_code)
       if (avgsize < mpio_config%output%par_min_block_size) then
          too_small=.true.
       else
          too_small=.false.
       end if
    end if
    call PAR_BROADCAST(too_small, wherein=wherein_code)
  end function is_too_small

  function is_parallel() result(paris_est_magique)
    logical(lg) paris_est_magique
    paris_est_magique=(mpio_config%output%post_process%parallel .and. kfl_reawr == 0) .or. (mpio_config%output%restart%parallel .and. (kfl_reawr == 2 .or. kfl_reawr == 1))
  end function is_parallel

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Destroy communicator
  !> @details Destroy postprocess communicator
  !> 
  !-----------------------------------------------------------------------

  subroutine posmpio_destructor()
    
    call comm_npoin_postpr_post % deallo()
    call comm_nelem_postpr_post % deallo()
    call comm_nboun_postpr_post % deallo()

    call comm_npoin_postpr_rst % deallo()
    call comm_nelem_postpr_rst % deallo()
    call comm_nboun_postpr_rst % deallo()

  end subroutine posmpio_destructor
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Redistribute
  !> @details Construct communicators for redistribution
  !> 
  !-----------------------------------------------------------------------

  subroutine posmpio_redistribute(FORCE_GENERATE_COMMUNICATOR)

    use def_domain
    use def_master
    use def_mpio
    use mod_redistribution, only : generate_comm_chunkdist
    use mod_redistribution, only : redistribution_array
    use mod_parall,         only : PAR_COMM_MY_CODE_WM
    use mod_communications, only : PAR_SUM
    use mod_communications, only : PAR_COMM_RANK_AND_SIZE
    !use mod_debug

    implicit none

    logical(lg), optional, intent(in) :: FORCE_GENERATE_COMMUNICATOR
    integer(ip)                       :: chunk_npoin_postpr
    integer(ip)                       :: chunk_nelem_postpr
    integer(ip)                       :: chunk_nboun_postpr
    integer(ip)                       :: chunk_npoin
    integer(ip)                       :: chunk_nelem
    integer(ip)                       :: chunk_nboun
    integer(ip)                       :: sum_npoin
    integer(ip)                       :: sum_nelem
    integer(ip)                       :: sum_nboun
    integer(4)                        :: iproc_par4,nproc_par4
    logical(lg)                       :: if_force_communicator

    if( INOTMASTER .and. mpio_config%output%merge ) then

       if( present(FORCE_GENERATE_COMMUNICATOR) ) then
          if_force_communicator = FORCE_GENERATE_COMMUNICATOR
       else
          if_force_communicator = .false.
       end if
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE_WM,iproc_par4,nproc_par4)
       !
       ! Evaluate chunks sizes
       !
       if( if_force_communicator .or. ((kfl_reawr==0 .and. .not.init_redistrib_post) .or. (kfl_reawr==2 .and. .not.init_redistrib_rst))) then
          
          sum_npoin  = meshe(imesh) % npoin_own
          sum_nelem  = meshe(imesh) % nelem
          sum_nboun  = meshe(imesh) % nboun

          call PAR_SUM(sum_npoin,PAR_COMM_MY_CODE_WM,INCLUDE_ROOT=.true.)
          call PAR_SUM(sum_nelem,PAR_COMM_MY_CODE_WM,INCLUDE_ROOT=.true.)
          call PAR_SUM(sum_nboun,PAR_COMM_MY_CODE_WM,INCLUDE_ROOT=.true.)
                       
          chunk_npoin_postpr = sum_npoin/int(nproc_par4,ip) + 1_ip
          chunk_nelem_postpr = sum_nelem/int(nproc_par4,ip) + 1_ip
          chunk_nboun_postpr = sum_nboun/int(nproc_par4,ip) + 1_ip

          chunk_npoin        = chunk_npoin_postpr
          chunk_nelem        = chunk_nelem_postpr
          chunk_nboun        = chunk_nboun_postpr

          if((iproc_par4+1)*chunk_npoin_postpr > sum_npoin ) then
             chunk_npoin=max(0_ip,sum_npoin-(iproc_par4)*chunk_npoin_postpr)
          end if
          if((iproc_par4+1)*chunk_nelem_postpr > sum_nelem ) then
             chunk_nelem=max(0_ip,sum_nelem-(iproc_par4)*chunk_nelem_postpr)
          end if
          if((iproc_par4+1)*chunk_nboun_postpr > sum_nboun ) then
             chunk_nboun=max(0_ip,sum_nboun-(iproc_par4)*chunk_nboun_postpr)
          end if
          !
          ! Evaluate my particular chunk size
          !
          if( kfl_reawr == 0 ) then

             init_redistrib_post = .true.
             chunk_npoin_post    = chunk_npoin
             chunk_nelem_post    = chunk_nelem
             chunk_nboun_post    = chunk_nboun

             call generate_comm_chunkdist(meshe(imesh) % lninv_loc,meshe(imesh) % npoin_own,chunk_npoin_postpr,PAR_COMM_MY_CODE_WM,&
                  &  comm_npoin_postpr_post,dist_dom='TARGET')
             call generate_comm_chunkdist(meshe(imesh) % leinv_loc,meshe(imesh) % nelem,chunk_nelem_postpr,PAR_COMM_MY_CODE_WM,&
                  &  comm_nelem_postpr_post,dist_dom='TARGET')
             call generate_comm_chunkdist(meshe(imesh) % lbinv_loc,meshe(imesh) % nboun,chunk_nboun_postpr,PAR_COMM_MY_CODE_WM,&
                  &  comm_nboun_postpr_post,dist_dom='TARGET')

          else if( kfl_reawr == 2 ) then

             init_redistrib_rst = .true.
             chunk_npoin_rst    = chunk_npoin
             chunk_nelem_rst    = chunk_nelem
             chunk_nboun_rst    = chunk_nboun

             call generate_comm_chunkdist(meshe(imesh) % lninv_loc,meshe(imesh) % npoin_own,chunk_npoin_postpr,PAR_COMM_MY_CODE_WM,&
                  &  comm_npoin_postpr_rst,dist_dom='TARGET')
             call generate_comm_chunkdist(meshe(imesh) % leinv_loc,meshe(imesh) % nelem,chunk_nelem_postpr,PAR_COMM_MY_CODE_WM,&
                  &  comm_nelem_postpr_rst,dist_dom='TARGET')
             call generate_comm_chunkdist(meshe(imesh) % lbinv_loc,meshe(imesh) % nboun,chunk_nboun_postpr,PAR_COMM_MY_CODE_WM,&
                  &  comm_nboun_postpr_rst,dist_dom='TARGET')
          end if

       end if

       if( kfl_reawr == 0 ) then

          comm_npoin_postpr => comm_npoin_postpr_post
          comm_nelem_postpr => comm_nelem_postpr_post
          comm_nboun_postpr => comm_nboun_postpr_post
          chunk_npoin       =  chunk_npoin_post
          chunk_nelem       =  chunk_nelem_post
          chunk_nboun       =  chunk_nboun_post

       else if( kfl_reawr == 2 ) then

          comm_npoin_postpr => comm_npoin_postpr_rst
          comm_nelem_postpr => comm_nelem_postpr_rst
          comm_nboun_postpr => comm_nboun_postpr_rst
          chunk_npoin       =  chunk_npoin_rst
          chunk_nelem       =  chunk_nelem_rst
          chunk_nboun       =  chunk_nboun_rst

       end if

       if( resultson(1:5) == 'NPOIN' ) then
          comm_postpr => comm_npoin_postpr
          dim         =  chunk_npoin
       else if( resultson(1:5) == 'NELEM' ) then
          comm_postpr => comm_nelem_postpr
          dim =          chunk_nelem
       else if( resultson(1:5) == 'NBOUN' ) then
          comm_postpr => comm_nboun_postpr
          dim         =  chunk_nboun
       end if

    end if

  end subroutine posmpio_redistribute

  subroutine redistribute_lnods()
    integer(ip)           :: ii,jj
    if (object(1:5)=='LNODS') then
       do ii=1, dim
          do jj=1, pdime_mpio
             if (ivec(jj,ii)>0) then
                ivec(jj,ii)=meshe(imesh) % lninv_loc(ivec(jj,ii))
             end if
          end do
       end do
    end if
  end subroutine redistribute_lnods

  subroutine redistribute_lelbo()
    integer(ip)           :: ii
    if (object(1:5)=='LELBO') then
       do ii=1, dim
          if (isca(ii)>0) then
             isca(ii)=meshe(imesh) % leinv_loc(isca(ii))
          end if
       end do
    end if
  end subroutine redistribute_lelbo

  subroutine redistribute_lnodb()
    integer(ip)           :: ii,jj
    if (object(1:5)=='LNODB') then
       do ii=1, dim
          do jj=1, pdime_mpio
             if (ivec(jj,ii)>0) then
                ivec(jj,ii)=meshe(imesh) % lninv_loc(ivec(jj,ii))
             end if
          end do
       end do
    end if
  end subroutine redistribute_lnodb

end module mod_mpio_par_postpr
