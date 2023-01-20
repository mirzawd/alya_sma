!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_partitioning.f90
!> @author  houzeaux
!> @date    2019-10-21
!> @brief   Parallel preprocess
!> @details Partitioning and redistribution to the slaves
!-----------------------------------------------------------------------

module mod_par_partitioning

  use def_kintyp,                      only : ip,lg,rp
  use def_master,                      only : IMASTER
  use def_master,                      only : kfl_paral
  use def_master,                      only : npart
  use def_master,                      only : routp
  use def_master,                      only : PART_AND_WRITE
  use def_master,                      only : READ_AND_RUN
  use def_master,                      only : PART_AND_RUN
  use mod_parall,                      only : PAR_PARALLEL_PARTITION
  use mod_parall,                      only : PAR_SEQUENTIAL_PARTITION
  use mod_parall,                      only : PAR_COMM_MY_CODE
  use def_parall,                      only : kfl_parseq_par 
  use def_parall,                      only : kfl_partition_par
  use def_parall,                      only : boxes_fine_par 
  use def_parall,                      only : boxes_coarse_par 
  use mod_parall,                      only : PAR_SFC
  use mod_parall,                      only : commd
  use mod_parall,                      only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                      only : PAR_COMM_MY_CODE
  use mod_parall,                      only : commd  
  use mod_par_parallel_partitioning,   only : par_parallel_partitioning
  use mod_par_sequential_partitioning, only : par_sequential_partitioning
  use mod_partition_sfc,               only : partition_sfc_statistics
  use mod_outfor,                      only : outfor
  use mod_par_parallel_restart,        only : par_parallel_restart
  use mod_par_additional_arrays,       only : par_ordered_exchange_update
  use mod_par_additional_arrays,       only : par_global_variables_arrays
  implicit none

  private

  public :: par_partitioning
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-10-21
  !> @brief   Partitioning
  !> @details Mesh partitioning
  !> 
  !-----------------------------------------------------------------------

  subroutine par_partitioning(READ_MESH,LCORRECT)
    
    logical(lg),       intent(in), optional  :: READ_MESH
    real(rp), pointer, intent(in), optional  :: LCORRECT(:)
    real(rp)                                 :: TIMINGS(5)

    !----------------------------------------------------------------------
    !
    ! Allocate communicator
    !
    !----------------------------------------------------------------------

    nullify(commd) 
    allocate(PAR_COMM_MY_CODE_ARRAY(1))
    call PAR_COMM_MY_CODE_ARRAY(1) % init(COMM_NAME='COMMD')
  
    PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD =  PAR_COMM_MY_CODE
    PAR_COMM_MY_CODE_ARRAY(1) % RANK4          =  int(kfl_paral,4)
    PAR_COMM_MY_CODE_ARRAY(1) % SIZE4          =  int(npart+1,4)
    commd                                      => PAR_COMM_MY_CODE_ARRAY(1)

    !----------------------------------------------------------------------
    !
    ! Partitioning
    !
    !----------------------------------------------------------------------

    select case ( kfl_parseq_par )

    case ( PAR_PARALLEL_PARTITION ) 
       !
       ! Parallel partitioning
       !
       if( READ_AND_RUN() ) then
          call par_parallel_restart()
       else
          call par_parallel_partitioning(READ_MESH,LCORRECT,TIMINGS=TIMINGS)
       end if

    case ( PAR_SEQUENTIAL_PARTITION ) 
       !
       ! Sequential partitioning
       !
       call par_sequential_partitioning()

    end select
    
    !----------------------------------------------------------------------
    !
    ! Output statistics
    !  
    !---------------------------------------------------------------------- 

    if( .not. READ_AND_RUN() ) then
       if( kfl_partition_par ==  PAR_SFC ) then
          call partition_sfc_statistics(&
               routp(1),routp(2),routp(3),routp(4),routp(5),boxes_fine_par,&
               boxes_coarse_par,COMM4=PAR_COMM_MY_CODE)
       end if
       call outfor(95_ip)
    end if

    call outfor(77_ip,REAL_NUMBERS=TIMINGS)

    !----------------------------------------------------------------------
    !
    ! Parallel preprocess, dump files
    !
    !----------------------------------------------------------------------

    if( kfl_parseq_par == PAR_PARALLEL_PARTITION .and. PART_AND_WRITE() ) then
       call par_parallel_restart()
       call runend('O.K.!')
    end if

    !----------------------------------------------------------------------
    !
    ! Others
    !
    !----------------------------------------------------------------------
    !
    ! Compute some useful variables
    !
    call par_global_variables_arrays()
    !
    ! Required dimensions
    !
    call par_mesh_dimensions()

  end subroutine par_partitioning

end module mod_par_partitioning
!> @}
