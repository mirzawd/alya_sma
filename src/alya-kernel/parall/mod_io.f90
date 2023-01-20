!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup IO
!> @{
!> @file    mod_io.f90
!> @author  houzeaux
!> @date    2020-02-24
!> @brief   Some IO tools
!> @details Tools for IO
!-----------------------------------------------------------------------

module mod_io 

  use def_kintyp,                    only : ip,rp,lg
  use def_master
  use def_domain,                    only : npoin
  use def_domain,                    only : nelem
  use def_domain,                    only : nboun
  use def_domain,                    only : memor_dom
  use mod_memory,                    only : memory_alloca
  use mod_memory,                    only : memory_deallo
  use mod_memory,                    only : memory_copy
  use mod_memory,                    only : memory_size
  use mod_redistribution,            only : redistribution_array
  use mod_communications,            only : PAR_SCATTER
  use mod_renumbering_nodes,         only : permr_nodes
  use mod_messages,                  only : messages_live
  use mod_parall,                    only : PAR_PARALLEL_PARTITION 
  use mod_parall,                    only : commd_npoin_from_mpio
  use mod_parall,                    only : commd_nelem_from_mpio
  use mod_parall,                    only : commd_nboun_from_mpio
  use def_parall,                    only : kfl_parseq_par
  use def_mpio,                      only : mpio_memor
  use mod_mpio_par_io,               only : PAR_FILE_READ_ALL
  use mod_mpio_par_readom,           only : npoin_vec
  use mod_mpio_par_readom,           only : nelem_vec
  use mod_mpio_par_readom,           only : nboun_vec
  use mod_communications_global,     only : PAR_MAX
  use mod_mpio_par_readom
  use mod_mpio_par_configure
  implicit none

  private

  interface io_read_merge
     module procedure &
          io_read_merge_rp_1,&
          io_read_merge_rp_2,&
          io_read_merge_ip_1         
  end interface io_read_merge
  
  integer(ip)               :: nenti_from_mpio 
  integer(ip),      pointer :: nenti_vec_ip(:)
  integer(ip)               :: nenti_size
  integer(ip)               :: nenti
  character(5)              :: wtype5
  character(150)            :: wherein
  
  public :: io_read_merge
   
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-24
  !> @brief   Compute dimensions
  !> @details Compute dimensions
  !> 
  !-----------------------------------------------------------------------
  
  subroutine io_dimensions(wtype)
    
    character(len=*), intent(in) :: wtype
    integer(ip)                 :: ii
    
    wherein = 'IN MPIO WITH MASTER'
    nullify(nenti_vec_ip)
    
    if( kfl_parseq_par /= PAR_PARALLEL_PARTITION ) call runend('YOU SHOULD USE SFC FOR READING PARALLEL FILES')
    !
    ! Check if *_VEC arrays exist
    !
    !if( .not. associated(npoin_vec) ) then
    !   call compute_dimensions()
    !   call par_init_share()
    !   call par_init_dimensions()
    !end if
    !
    ! Define sizes
    !
    wtype5  = wtype(1:5)
    if(      wtype5 == 'NPOIN' ) then
       nenti      = npoin
       nenti_size = size(npoin_vec,KIND=ip)
       if( .not. associated(npoin_vec)) call runend('IO_READ_MERGE: NPOIN_VEC IS NOT ASSOCIATED')
    else if( wtype5 == 'NELEM' ) then
       nenti      = nelem
       nenti_size = size(nelem_vec,KIND=ip)
       if( .not. associated(nelem_vec)) call runend('IO_READ_MERGE: NPOIN_VEC IS NOT ASSOCIATED')
    else if( wtype5 == 'NBOUN' ) then
       nenti      = nboun
       nenti_size = size(nboun_vec,KIND=ip)
       if( .not. associated(nboun_vec)) call runend('IO_READ_MERGE: NPOIN_VEC IS NOT ASSOCIATED')
    end if
    !
    ! Scatter dimensions
    ! copy the dimensions into the array of precision IP because npoin_vec is int(4) PAR_SCATTER expects int(ip)
    !    
    call memory_alloca(memor_dom,'NENTI_VEC_IP','io_read_merge',nenti_vec_ip,nenti_size)
    if(      wtype5 == 'NPOIN' ) then
       do ii = 1,nenti_size
          nenti_vec_ip(ii) = int(npoin_vec(ii),ip)
       end do
    else if( wtype5 == 'NELEM' ) then
       do ii = 1,nenti_size
          nenti_vec_ip(ii) = int(nelem_vec(ii),ip)
       end do
    else if( wtype5 == 'NBOUN' ) then
       do ii = 1,nenti_size
          nenti_vec_ip(ii) = int(nboun_vec(ii),ip)
       end do
    end if
    call PAR_SCATTER(nenti_vec_ip,nenti_from_mpio,wherein='IN MY CODE')
 
  end subroutine io_dimensions

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-24
  !> @brief   Read and redistribute
  !> @details Read and redistribute twice
  !> 
  !-----------------------------------------------------------------------

  subroutine io_read_merge_rp_1(xx,file_name,wtype)

    real(rp),         pointer, intent(inout) :: xx(:)
    character(len=*),          intent(in)    :: file_name
    character(len=*),          intent(in)    :: wtype
    integer(ip)                              :: ii
    real(rp),         pointer                :: xx_in(:)
    
    nullify(xx_in)
    !
    ! Number of entities
    !
    call io_dimensions(wtype)
    !
    ! Read with MPIO
    !
    call memory_alloca(memor_dom,'XX_IN','io_read_merge',xx_in,max(nenti_from_mpio,1_ip))
    call PAR_FILE_READ_ALL(xx_in, file_name, nenti_from_mpio, wherein=wherein)
    !
    ! Redistribute after reading
    !
    if( INOTMASTER ) then
       if(      wtype5 == 'NPOIN' ) then
          call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',COMM=commd_npoin_from_mpio,PERMUTE_RECV=.true.,VERBOSE=.false.)
       else if( wtype5 == 'NELEM' ) then
          continue
       else if( wtype5 == 'NBOUN' ) then
          call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',COMM=commd_nboun_from_mpio,PERMUTE_RECV=.false.,VERBOSE=.false.)
       end if
    end if
    !
    ! Redistribute using partitioning communicator
    !
    call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',VERBOSE=.false.)
    !
    ! Copy
    !
    do ii = 1,nenti
       xx(ii) = xx_in(ii)
    end do
    
    call memory_deallo(memor_dom,'XX_IN'       ,'io_read_merge',xx_in)
    call memory_deallo(memor_dom,'NENTI_VEC_IP','io_read_merge',nenti_vec_ip)

  end subroutine io_read_merge_rp_1

  subroutine io_read_merge_rp_2(xx,file_name,wtype)

    real(rp),         pointer, intent(inout) :: xx(:,:)
    character(len=*),          intent(in)    :: file_name
    character(len=*),          intent(in)    :: wtype
    integer(ip)                              :: ii,kdime
    real(rp),         pointer                :: xx_in(:,:)
    
    nullify(xx_in)
    kdime = memory_size(xx,1_ip)
    call PAR_MAX(kdime)
    !
    ! Number of entities
    !
    call io_dimensions(wtype)
    !
    ! Read with MPIO
    !
    call memory_alloca(memor_dom,'XX_IN','io_read_merge',xx_in,kdime,max(nenti_from_mpio,1_ip))
    call PAR_FILE_READ_ALL(xx_in, file_name, kdime, nenti_from_mpio, wherein=wherein)
    !
    ! Redistribute after reading
    !
    if( INOTMASTER ) then       
       if(      wtype5 == 'NPOIN' ) then
          call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',COMM=commd_npoin_from_mpio,PERMUTE_RECV=.true.,VERBOSE=.false.)
       else if( wtype5 == 'NELEM' ) then
          continue
       else if( wtype5 == 'NBOUN' ) then
          call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',COMM=commd_nboun_from_mpio,PERMUTE_RECV=.false.,VERBOSE=.false.)
       end if
    end if
    !
    ! Redistribute using partitioning communicator
    !
    call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',VERBOSE=.false.)
    !
    ! Copy
    !
    do ii = 1,nenti
       xx(1:kdime,ii) = xx_in(1:kdime,ii)
    end do
    
    call memory_deallo(memor_dom,'XX_IN'       ,'io_read_merge',xx_in)
    call memory_deallo(memor_dom,'NENTI_VEC_IP','io_read_merge',nenti_vec_ip)

  end subroutine io_read_merge_rp_2

    subroutine io_read_merge_ip_1(xx,file_name,wtype)

    integer(ip),      pointer, intent(inout) :: xx(:)
    character(len=*),          intent(in)    :: file_name
    character(len=*),          intent(in)    :: wtype
    integer(ip)                              :: ii
    integer(ip),      pointer                :: xx_in(:)
    
    nullify(xx_in)
    !
    ! Number of entities
    !
    call io_dimensions(wtype)
    !
    ! Read with MPIO
    !
    call memory_alloca(memor_dom,'XX_IN','io_read_merge',xx_in,max(nenti_from_mpio,1_ip))
    call PAR_FILE_READ_ALL(xx_in, file_name, nenti_from_mpio, wherein=wherein)
    !
    ! Redistribute after reading
    !
    if( INOTMASTER ) then
       if(      wtype5 == 'NPOIN' ) then
          call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',COMM=commd_npoin_from_mpio,PERMUTE_RECV=.true.,VERBOSE=.false.)
       else if( wtype5 == 'NELEM' ) then
          continue
       else if( wtype5 == 'NBOUN' ) then
          call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',COMM=commd_nboun_from_mpio,PERMUTE_RECV=.false.,VERBOSE=.false.)
       end if
    end if
    !
    ! Redistribute using partitioning communicator
    !
    call redistribution_array(xx_in,wtype5,MEMOR=memor_dom,VARIABLE_NAME='XX_IN',VERBOSE=.false.)
    !
    ! Copy
    !
    do ii = 1,nenti
       xx(ii) = xx_in(ii)
    end do
    
    call memory_deallo(memor_dom,'XX_IN'       ,'io_read_merge',xx_in)
    call memory_deallo(memor_dom,'NENTI_VEC_IP','io_read_merge',nenti_vec_ip)

  end subroutine io_read_merge_ip_1

end module mod_io
!> @}
