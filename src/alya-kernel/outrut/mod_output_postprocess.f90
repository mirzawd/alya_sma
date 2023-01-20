!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    mod_postprocess.f90
!> @author  houzeaux
!> @date    2018-11-08
!> @brief   Postprocess sruff
!> @details Output and postprocess stuffs
!-----------------------------------------------------------------------

module mod_output_postprocess

  use def_parame
  use def_domain
  use def_master 
  use def_inpout
  use def_postpr
  use def_kermod,         only : nwith
  use def_kermod,         only : witness_mesh
  use mod_memory,         only : memory_alloca
  use mod_communications, only : PAR_SUM
  use mod_ecoute,         only : ecoute
  use mod_witness,        only : WITNESS_INSTANTANEAOUS
  use mod_witness,        only : WITNESS_AVERAGE
  use mod_witness,        only : WITNESS_ACCUMULATE
  use mod_witness,        only : WITNESS_SPHERE
  use mod_witness,        only : WITNESS_BOX
  use mod_witness,        only : WITNESS_RING
  use mod_witness,        only : WITNESS_BOUNDARY
  use mod_witness,        only : WITNESS_ELEMENT_SET
  use mod_witness,        only : WITNESS_MATERIAL 
  use mod_witness,        only : WITNESS_PLANE 
  use mod_witness,        only : WITNESS_BOUNDARY_SET
  use mod_witness,        only : WITNESS_BAR3D
  use mod_witness,        only : WITNESS_3D_CIRCLE
  use mod_witness,        only : WITNESS_3D_RING
  use mod_witness,        only : WITNESS_IN_BOX
  use mod_witness,        only : witness_geometry_averaging_ini
  use mod_witness,        only : witness_geometry_averaging_end
  use mod_witness,        only : witness_point_averaging_ini
  use mod_witness,        only : witness_point_averaging_end
  use mod_exchange,       only : exchange_add
  use mod_maths_basic,    only : maths_normalize_vector
  
  implicit none
  
  abstract interface
     subroutine xxx_outvar(ivari,imesh)
#ifdef I8
       integer(8) :: ivari
       integer(8) :: imesh
#else
       integer(4) :: ivari
       integer(4) :: imesh
#endif
     end subroutine xxx_outvar
  end interface

  character(5) :: CNULL='NULL '

  private

  public :: output_postprocess_initialization                 
  public :: output_postprocess_allocate                       
  public :: output_postprocess_parall                         
  public :: output_postprocess_parall_old                     
  public :: output_postprocess_read                           
  public :: output_postprocess_allocate_sets_and_witness      
  public :: output_postprocess_header_sets_and_witness   
  public :: output_postprocess_output_witness                 
  public :: output_postprocess_element_sets_parall            
  public :: output_postprocess_boundary_sets_parall           
  public :: output_postprocess_node_sets_parall               
  public :: output_postprocess_witness_parall                 
  public :: output_postprocess_cancel_variable_postprocess    
  public :: output_postprocess_sets_ordering                  
  public :: output_postprocess_check_witness_output           
  public :: output_postprocess_check_variable_postprocess_now 
  public :: output_postprocess_check_variable_postprocess     
  public :: output_postprocess_check_variable_witness
  public :: output_postprocess_check_variable_node_sets
  public :: output_postprocess_check_variable_boundary_sets
  public :: output_postprocess_check_variable_element_sets
  public :: output_postprocess_read_geometry                  ! Read a geometry
  public :: output_postprocess_variables                      ! Postprocess variables of a module
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Initialization
  !> @details Allocation and inititalization of POSTP structure
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_initialization()

    integer(ip) :: imodu
    
    do imodu = 0,mmodu
       nullify( momod(imodu) % postp )
    end do

  end subroutine output_postprocess_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Allocate
  !> @details Allocation and inititalization of POSTP structure
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_allocate(imodu)

    integer(ip), optional :: imodu
    integer(ip)           :: kmodu,ivarp,ivars,ivarw,ivari
    integer(ip)           :: imodu_ini,imodu_end

    if( present(imodu) ) then
       imodu_ini = imodu
       imodu_end = imodu
    else
       imodu_ini = 1
       imodu_end = mmodu
    end if

    do kmodu = imodu_ini,imodu_end

       if( kfl_modul(kmodu) /= 0 .and. ( .not. associated(momod(kmodu) % postp) ) )  then
          !
          ! Allocate 
          !
          allocate( momod(kmodu) % postp(1) )
          !
          ! Postprocess on mesh
          !
          postp => momod(kmodu)%postp
          postp(1) % npp_inits  = 0                          ! Postprocess initial step
          postp(1) % kfl_oonce  = 0                          ! Do not postprocess only once
          postp(1) % npp_stepi  = 0                          ! Postprocess step interval
          postp(1) % pos_alrea  = 0                          ! Already postprocessed
          postp(1) % vox_stepi  = 0                          ! Postprocess step interval
          postp(1) % vox_alrea  = 0                          ! Already postprocessed
          postp(1) % npp_iniso  = 0                          ! Postprocess initial condition
          postp(1) % pos_tinit  = 0.0_rp                     ! Postprocess initial time
          postp(1) % pos_times  = 0.0_rp                     ! Postprocess times 
          postp(1) % vox_times  = 0.0_rp                     ! Voxel times 
          postp(1) % pos_perio  = 0.0_rp                     ! Postprocess time period
          postp(1) % vox_perio  = 0.0_rp                     ! Voxel time period
          postp(1) % lcomp      = -1                         ! Postprocess all components
          postp(1) % rst_time   = .true.
          !
          ! Postprocess on sets, wtiness and co
          !
          postp(1) % npp_setse     = 0                          ! Postprocess element sets calculation
          postp(1) % npp_setsb     = 0                          ! Postprocess boundary sets calculation
          postp(1) % npp_setsn     = 0                          ! Postprocess node sets calculation
          postp(1) % per_setse     = 0                          ! Postprocess permutation element sets calculation
          postp(1) % per_setsb     = 0                          ! Postprocess permutation boundary sets calculation
          postp(1) % per_setsn     = 0                          ! Postprocess permutation node sets calculation
          postp(1) % npp_stepelset = 1                          ! Postprocess element sets at every step
          postp(1) % npp_stepnoset = 1                          ! Postprocess node sets at every step
          postp(1) % npp_stepboset = 1                          ! Postprocess boundary sets at every step
          postp(1) % npp_stepw     = 1                          ! Postprocess witness points at every step
          postp(1) % npp_stepg     = 0                          ! Postprocess witness points at every step
          postp(1) % npp_witne     = 0                          ! Postprocess witness points
          postp(1) % npp_witng     = 0                          ! Postprocess witness geometries
          do ivarw = 1,nvarg
             postp(1) % witng_kfl_time(ivarw)  = WITNESS_INSTANTANEAOUS  ! Witness geometry: Time strategy
             postp(1) % witng_perio_max(ivarw) = 0.0_rp                  ! Witness geometry: Averaging period
          end do
          do ivarw = 1,nvarw
             postp(1) % witne_kfl_time(ivarw)= WITNESS_INSTANTANEAOUS  ! Witness point: Time strategy
             postp(1) % witne_period(ivarw)  = 0.0_rp                  ! Witness point: Averaging period
             postp(1) % witne_average(ivarw) = 0.0_rp                  ! Witness point: Averaging time
             postp(1) % witne_dt(ivarw)      = 1.0_rp                  ! Witness point: Averaging dt
          end do

          postp(1) % nvaes     = 0                           ! Element  set variables
          postp(1) % nvabs     = 0                           ! Boundary set variables
          postp(1) % nvans     = 0                           ! Node     set variables
          postp(1) % per_nvaes = 0                           ! Element  set variables
          postp(1) % per_nvabs = 0                           ! Boundary set variables
          postp(1) % per_nvans = 0                           ! Node     set variables
          postp(1) % nvawi     = 0                           ! Witness point variables
          postp(1) % nvawg     = 0                           ! Witness geometry variables
          postp(1) % lun_setse = 0                           ! Element set unit imodu*10+6
          postp(1) % lun_setsb = 0                           ! Boundary set unit imodu*10+7
          postp(1) % lun_setsn = 0                           ! Node set unit imodu*10+8
          postp(1) % lun_witne = 0                           ! Node set unit imodu*10+8
          postp(1) % ipass     = 0                           ! Set memory allocated and header
          !
          ! Real(rp)
          !
          do ivars = 1,nvars
             do ivari = 1,5
                postp(1) % paese(ivari,ivars) = 0.0_rp       ! Element set parameters  
                postp(1) % pabse(ivari,ivars) = 0.0_rp       ! Boundary sets parameters  
                postp(1) % panse(ivari,ivars) = 0.0_rp       ! Node set parameters  
             end do
          end do
          !
          ! Character(5)
          !
          do ivars = 1,nvars
             postp(1) % woese(ivars) = cnull                 ! Name of the element set variables
             postp(1) % wobse(ivars) = cnull                 ! Name of the boundary set variables
             postp(1) % wonse(ivars) = cnull                 ! Name of the node set variables          
          end do
          do ivarw = 1,nvarw
             postp(1) % wowit(ivarw) = cnull                 ! Name of the witness point variables
          end do
          do ivarw = 1,nvarg
             postp(1) % wowig(ivarw) = cnull                 ! Name of the witness geometry variables
          end do
          do ivarp = 1,nvarp
             postp(1) % wopos(1,ivarp) = cnull               ! Name of postprocess variable
             postp(1) % wopos(2,ivarp) = 'SCALA'             ! Character of postprocess variable
             postp(1) % wopos(3,ivarp) = 'NPOIN'             ! NPOIN type by default
          end do
          !
          ! Register
          !
          do ivarp = 1,nvarp
             postp(1) % enti_posit(ivarp)       = 0
             postp(1) % comp_posit(ivarp)       = 0
             postp(1) % time_posit(ivarp)       = 0
             postp(1) % time_num(ivarp)         = 0
             postp(1) % comp_num(ivarp)         = 0
             postp(1) % dime_num(ivarp)         = 0         
             postp(1) % array_allocated(ivarp)  = 0
             postp(1) % array_used(ivarp)       = 0
             postp(1) % array_registered(ivarp) = 0
          end do
          !
          ! Pointers
          !
          nullify( postp(1) % veset )
          nullify( postp(1) % vbset )
          nullify( postp(1) % vnset )
          nullify( postp(1) % witne )
          nullify( postp(1) % witng )
          nullify( postp(1) % witng_perio_count )
          nullify( postp(1) % witng_denom )
          nullify( postp(1) % witng_deldenom )
          nullify( postp(1) % witng_delperio )

       end if

    end do

  end subroutine output_postprocess_allocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Parallelization
  !> @details Exchange POSTP structure 
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_parall_old()

    integer(ip) :: ivarp,ivars,ivarw,ivart,ivari,imesh

    postp => momod(modul)%postp
    !
    ! Mesh postprocess
    !
    postp => momod(modul)%postp
    call iexcha(postp(1) % npp_inits)                             ! Postprocess initial step
    call iexcha(postp(1) % npp_iniso)                             ! Postprocess initial condition
    call rexcha(postp(1) % pos_tinit)                             ! Postprocess initial time
    do imesh = 0,nvarm
       do ivarp = 1,nvarp
          call iexcha(postp(1) % npp_stepi(ivarp,imesh))          ! Postprocess step interval
          call iexcha(postp(1) % pos_alrea(ivarp,imesh))          ! Already postprocessed
          call iexcha(postp(1) % vox_stepi(ivarp,imesh))          ! Postprocess step interval
          call iexcha(postp(1) % vox_alrea(ivarp,imesh))          ! Already postprocessed
          call rexcha(postp(1) % pos_perio(ivarp,imesh))          ! Postprocess time period 
          call rexcha(postp(1) % vox_perio(ivarp,imesh))          ! Postprocess time period 
       end do
    end do
    do imesh = 0,nvarm
       do ivarp = 1,nvarp
          do ivari = 1,ncomp_max
             call iexcha(postp(1) % lcomp(ivari,ivarp,imesh))     ! List of components to postprocess
          end do
       end do
    end do
    do ivarp = 1,nvarp
       call iexcha(postp(1) % kfl_oonce(ivarp))                   ! Postprocess only once
    end do

    do imesh = 0,nvarm
       do ivarp = 1,nvarp
          do ivart = 1,nvart
             call rexcha(postp(1) % pos_times(ivart,ivarp,imesh)) ! Postprocess times 
             call rexcha(postp(1) % vox_times(ivart,ivarp,imesh)) ! Postprocess times 
          end do
       end do
    end do
    !
    ! Array definition
    !
    do ivarp = 1,nvarp
       call iexcha(postp(1) % enti_posit(ivarp))                ! Entity position
       call iexcha(postp(1) % comp_posit(ivarp))                ! Component position
       call iexcha(postp(1) % time_posit(ivarp))                ! Time position        
    end do
    !
    ! Sets, wtiness and co
    !
    do ivars = 1,nvars
       call iexcha(postp(1) % npp_setse(ivars))                 ! Postprocess element sets calculation
       call iexcha(postp(1) % npp_setsb(ivars))                 ! Postprocess boundary sets calculation
       call iexcha(postp(1) % npp_setsn(ivars))                 ! Postprocess node sets calculation
       call iexcha(postp(1) % per_setse(ivars))                 ! Postprocess element sets calculation
       call iexcha(postp(1) % per_setsb(ivars))                 ! Postprocess boundary sets calculation
       call iexcha(postp(1) % per_setsn(ivars))                 ! Postprocess node sets calculation
    end do
    call iexcha(postp(1) % npp_stepelset)                       ! Postprocess element sets interval
    call iexcha(postp(1) % npp_stepnoset)                       ! Postprocess node sets interval
    call iexcha(postp(1) % npp_stepboset)                       ! Postprocess boundary sets interval
    call iexcha(postp(1) % npp_stepw)                           ! Postprocess witness points interval
    call iexcha(postp(1) % npp_stepg)                           ! Postprocess witness geometries interval
    do ivarw = 1,nvarw
       call iexcha(postp(1) % npp_witne(ivarw))                 ! Postprocess witness points
    end do
    do ivarw = 1,nvarg
       call iexcha(postp(1) % npp_witng(ivarw))                 ! Postprocess witness geometries
    end do
    do ivarw = 1,nvarw
       call iexcha(postp(1) % witne_kfl_time(ivarw))            ! Witness point: Time strategy
       call rexcha(postp(1) % witne_period(ivarw))              ! Witness point: Averaging period
       call rexcha(postp(1) % witne_average(ivarw))             ! Witness point: Averaging time
       call rexcha(postp(1) % witne_dt(ivarw))                  ! Witness point: Time step
    end do
    do ivarw = 1,nvarg
       call iexcha(postp(1) % witng_kfl_time(ivarw))            ! Witness geometry: Time strategy
       call rexcha(postp(1) % witng_perio_max(ivarw))           ! Witness geometry: Averaging period
    end do

    call iexcha(postp(1) % nvaes)                               ! Element  set variables  
    call iexcha(postp(1) % nvabs)                               ! Boundary set variables
    call iexcha(postp(1) % nvans)                               ! Node set variables
    call iexcha(postp(1) % per_nvaes)                           ! Element  set variables  
    call iexcha(postp(1) % per_nvabs)                           ! Boundary set variables
    call iexcha(postp(1) % per_nvans)                           ! Node set variables
    call iexcha(postp(1) % nvawi)                               ! Witness point variables
    call iexcha(postp(1) % nvawg)                               ! Witness geometry variables
    call iexcha(postp(1) % lun_setse)                           ! Element set unit imodu*10+6
    call iexcha(postp(1) % lun_setsb)                           ! Boundary set unit imodu*10+7
    call iexcha(postp(1) % lun_setsn)                           ! Node set unit imodu*10+8
    call iexcha(postp(1) % lun_witne)                           ! Node set unit imodu*10+8
    call iexcha(postp(1) % lun_witng)                           ! Node set unit imodu*10+9
    call iexcha(postp(1) % ipass)                               ! Set memory allocated and header
    !
    ! Real(rp)
    !
    do ivars = 1,nvars
       do ivari = 1,5
          call rexcha(postp(1) % paese(ivari,ivars))            ! Element set parameters  
          call rexcha(postp(1) % pabse(ivari,ivars))            ! Boundary sets parameters  
          call rexcha(postp(1) % panse(ivari,ivars))            ! Node set parameters  
       end do
    end do
    !
    ! Character(5)
    !
    do ivars = 1,nvars
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % woese(ivars)
       if(parii==2.and.ISLAVE)  postp(1) % woese(ivars) = parch(nparc+1:nparc+5)
       nparc = nparc+5
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % wobse(ivars)
       if(parii==2.and.ISLAVE)  postp(1) % wobse(ivars) = parch(nparc+1:nparc+5)
       nparc = nparc+5
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % wonse(ivars)
       if(parii==2.and.ISLAVE)  postp(1) % wonse(ivars) = parch(nparc+1:nparc+5)
       nparc = nparc+5
    end do
    do ivarw = 1,nvarw
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % wowit(ivarw)
       if(parii==2.and.ISLAVE)  postp(1) % wowit(ivarw) = parch(nparc+1:nparc+5)
       nparc = nparc+5           
    end do
    do ivarp = 1,nvarp
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(1,ivarp) 
       if(parii==2.and.ISLAVE)  postp(1) % wopos(1,ivarp) = parch(nparc+1:nparc+5)
       nparc = nparc+5           
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(2,ivarp) 
       if(parii==2.and.ISLAVE)  postp(1) % wopos(2,ivarp) = parch(nparc+1:nparc+5)
       nparc = nparc+5           
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(3,ivarp) 
       if(parii==2.and.ISLAVE)  postp(1) % wopos(3,ivarp) = parch(nparc+1:nparc+5)
       nparc = nparc+5           
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(4,ivarp) 
       if(parii==2.and.ISLAVE)  postp(1) % wopos(4,ivarp) = parch(nparc+1:nparc+5)
       nparc = nparc+5           
       if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(5,ivarp) 
       if(parii==2.and.ISLAVE)  postp(1) % wopos(5,ivarp) = parch(nparc+1:nparc+5)
       nparc = nparc+5           
    end do

    if( nparc > len(parch) ) call runend('MOD_OUTPUT_POTPROCESS: TOO MANY CHARACTERS')

  end subroutine output_postprocess_parall_old

  subroutine output_postprocess_parall(postp_in)

    type(typos), optional, pointer, intent(inout) :: postp_in(:)
    type(typos),           pointer                :: postp_loc(:)
    integer(ip)                                   :: ivarp,ivars,ivarw

    if( present(postp_in) ) then
       postp_loc => postp_in
    else
       postp_loc => momod(modul) % postp
    end if
    !
    ! Mesh postprocess
    !
    call exchange_add(postp_loc(1) % npp_inits)                           ! Postprocess initial step
    call exchange_add(postp_loc(1) % npp_iniso)                           ! Postprocess initial condition
    call exchange_add(postp_loc(1) % pos_tinit)                           ! Postprocess initial time
    call exchange_add(postp_loc(1) % npp_stepi)                           ! Postprocess step interval
    call exchange_add(postp_loc(1) % pos_alrea)                           ! Already postp_locrocessed
    call exchange_add(postp_loc(1) % vox_stepi)                           ! Postprocess step interval
    call exchange_add(postp_loc(1) % vox_alrea)                           ! Already postp_locrocessed
    call exchange_add(postp_loc(1) % pos_perio)                           ! Postprocess time period 
    call exchange_add(postp_loc(1) % vox_perio)                           ! Postprocess time period 
    call exchange_add(postp_loc(1) % lcomp)                               ! List of components to postp_locrocess
    call exchange_add(postp_loc(1) % kfl_oonce)                           ! Postprocess only once

    call exchange_add(postp_loc(1) % pos_times)                           ! Postprocess times 
    call exchange_add(postp_loc(1) % vox_times)                           ! Postprocess times 
    !
    ! Array definition
    !
    call exchange_add(postp_loc(1) % enti_posit)                          ! Entity position
    call exchange_add(postp_loc(1) % comp_posit)                          ! Component position
    call exchange_add(postp_loc(1) % time_posit)                          ! Time position        
    !
    ! Sets, wtiness and co
    !
    call exchange_add(postp_loc(1) % npp_setse)                           ! Postprocess element sets calculation
    call exchange_add(postp_loc(1) % npp_setsb)                           ! Postprocess boundary sets calculation
    call exchange_add(postp_loc(1) % npp_setsn)                           ! Postprocess node sets calculation
    call exchange_add(postp_loc(1) % per_setse)                           ! Postprocess element sets calculation
    call exchange_add(postp_loc(1) % per_setsb)                           ! Postprocess boundary sets calculation
    call exchange_add(postp_loc(1) % per_setsn)                           ! Postprocess node sets calculation
    call exchange_add(postp_loc(1) % npp_stepelset)                       ! Postprocess element set interval
    call exchange_add(postp_loc(1) % npp_stepnoset)                       ! Postprocess element set interval
    call exchange_add(postp_loc(1) % npp_stepboset)                       ! Postprocess element set interval
    call exchange_add(postp_loc(1) % npp_stepw)                           ! Postprocess witness points interval
    call exchange_add(postp_loc(1) % npp_stepg)                           ! Postprocess witness geometries interval
    call exchange_add(postp_loc(1) % npp_witne)                           ! Postprocess witness points
 
    call exchange_add(postp_loc(1) % npp_witng)                           ! Postprocess witness geometries

    call exchange_add(postp_loc(1) % witne_kfl_time)                      ! Witness point: Time strategy
    call exchange_add(postp_loc(1) % witne_period)                        ! Witness point: Averaging period
    call exchange_add(postp_loc(1) % witne_average)                       ! Witness point: Averaging time
    call exchange_add(postp_loc(1) % witne_dt)                            ! Witness point: Time step
    
    call exchange_add(postp_loc(1) % witng_kfl_time)                      ! Witness geometry: Time strategy
    call exchange_add(postp_loc(1) % witng_perio_max)                     ! Witness geometry: Averaging period

    call exchange_add(postp_loc(1) % nvaes)                               ! Element  set variables  
    call exchange_add(postp_loc(1) % nvabs)                               ! Boundary set variables
    call exchange_add(postp_loc(1) % nvans)                               ! Node set variables
    call exchange_add(postp_loc(1) % per_nvaes)                           ! Element  set variables  
    call exchange_add(postp_loc(1) % per_nvabs)                           ! Boundary set variables
    call exchange_add(postp_loc(1) % per_nvans)                           ! Node set variables
    call exchange_add(postp_loc(1) % nvawi)                               ! Witness point variables
    call exchange_add(postp_loc(1) % nvawg)                               ! Witness geometry variables
    call exchange_add(postp_loc(1) % lun_setse)                           ! Element set unit imodu*10+6
    call exchange_add(postp_loc(1) % lun_setsb)                           ! Boundary set unit imodu*10+7
    call exchange_add(postp_loc(1) % lun_setsn)                           ! Node set unit imodu*10+8
    call exchange_add(postp_loc(1) % lun_witne)                           ! Node set unit imodu*10+8
    call exchange_add(postp_loc(1) % lun_witng)                           ! Node set unit imodu*10+9
    call exchange_add(postp_loc(1) % ipass)                               ! Set memory allocated and header
    !
    ! Real(rp)
    !
    call exchange_add(postp_loc(1) % paese)                               ! Element set parameters  
    call exchange_add(postp_loc(1) % pabse)                               ! Boundary sets parameters  
    call exchange_add(postp_loc(1) % panse)                               ! Node set parameters  
    !
    ! Character(5)
    !
    do ivars = 1,nvars
       call exchange_add(postp_loc(1) % woese(ivars))
       call exchange_add(postp_loc(1) % wobse(ivars))
       call exchange_add(postp_loc(1) % wonse(ivars))
    end do
    do ivarw = 1,nvarw       
       call exchange_add(postp_loc(1) % wowit(ivarw))
    end do
    do ivarp = 1,nvarp    
       call exchange_add(postp_loc(1) % wopos(1,ivarp))
       call exchange_add(postp_loc(1) % wopos(2,ivarp))
       call exchange_add(postp_loc(1) % wopos(3,ivarp))
       call exchange_add(postp_loc(1) % wopos(4,ivarp))
       call exchange_add(postp_loc(1) % wopos(5,ivarp))
    end do
    
  end subroutine output_postprocess_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Read
  !> @details Read from data file the postprocess and output options
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_read()

    integer(ip)                :: ivarp,ivarw,ivart,ivars
    integer(ip)                :: ipost,nvabs,nvaes,ii,iposi,imesh
    character(5)               :: mesh_name
    
    if( words(1) == 'POSTP' ) then
       if( exists('ONVOX') ) then
          ipost = 2
       else
          ipost = 1
       end if
       
       if( words(2) == 'INITI' ) then
          !
          ! Initial solution
          !
          postp(1) % npp_iniso = 1
          
       else
          
          do ivarp = 1,nvarp
             if( words(2) == postp(1) % wopos(1,ivarp) ) then
                iposi = 2
                !
                ! Mesh
                !                
                if( exists('MESH ') ) then
                   mesh_name = getcha('MESH ','NULL ','#Postprocess mesh')
                   if( associated(witness_mesh) ) then
                      do ii = 1,size(witness_mesh)
                         if( trim(mesh_name) == witness_mesh(ii) % name(1:5) ) then
                            imesh = ii
                         end if
                      end do
                   else
                      call runend('MOD_OUTPUT_POSTPROCESS: DEFINE WITNESS MESH '//trim(mesh_name))
                   end if
                else
                   imesh = 0
                end if
                !
                ! Components
                !                
                if( exists('COMPO') ) then
                   if( getcha('COMPO','ALL  ','#List of components to postprocess') == 'ALL  ') then
                      postp(1) % lcomp (:,ivarp,imesh) = -1_ip
                   else
                      ii = 1
                      do while( int(param(iposi+ii),ip) /= 0 )
                         postp(1) % lcomp (ii,ivarp,imesh) =  int(param(iposi+ii),ip)
                         ii = ii + 1
                      end do
                   end if
                end if
                
                if( exists('ONLYO') ) then
                   postp(1) % kfl_oonce(ivarp) = 1
                end if
                if( exists('NOTON') ) then
                   postp(1) % kfl_oonce(ivarp) = 0
                end if

                if( words(3) == 'STEPS' ) then
                   !
                   ! Steps
                   !
                   if( ipost == 1 ) then
                      postp(1) % npp_stepi(ivarp,imesh) = &
                           getint('STEPS',1_ip,&
                           '#Postprocess step interval for '// postp(1) % wopos(1,ivarp))
                      if (npp_stepo>-1) postp(1) % npp_stepi(ivarp,imesh) = npp_stepo
                   else
                      postp(1) % vox_stepi(ivarp,imesh) = &
                           getint('STEPS',1_ip,&
                           '#Postprocess step interval for '// postp(1) % wopos(1,ivarp))                       
                      if (npp_stepo>-1) postp(1) % vox_stepi(ivarp,imesh) = npp_stepo
                   end if
                   iposi = iposi + 1
                   if(exists('ATTIM')) then
                      !
                      ! Specific times
                      !
                      if( ipost == 1 ) then
                         do ivart = 1,nvart
                            postp(1) % pos_times(ivart,ivarp,imesh) = param(ivart+iposi)
                         end do
                      else
                         do ivart = 1,nvart
                            postp(1) % vox_times(ivart,ivarp,imesh) = param(ivart+iposi)
                         end do
                      end if
                      iposi = iposi + nvart
                   end if
                else if(words(3)=='PERIO') then
                   !
                   ! Time period
                   !
                   if( ipost == 1 ) then
                      postp(1) % pos_perio(ivarp,imesh)= getrea('PERIO',0.0_rp,'#Postprocess time period')
                   else
                      postp(1) % vox_perio(ivarp,imesh)= getrea('PERIO',0.0_rp,'#Postprocess time period')
                   end if
                   iposi = iposi + 1
                else
                   if( ipost == 1 ) then
                      postp(1) % npp_stepi(ivarp,imesh) = 1
                      if (npp_stepo>-1) postp(1) % npp_stepi(ivarp,imesh) = npp_stepo
                   else
                      postp(1) % vox_stepi(ivarp,imesh) = 1
                      if (npp_stepo>-1) postp(1) % vox_stepi(ivarp,imesh) = npp_stepo
                   end if
                end if
             end if
          end do
       end if

    else if( words(1) == 'ELEME' ) then
       !
       ! Element sets
       !
       postp(1) % npp_stepelset = getint('STEPS',1_ip,'#When to postprocess element sets')

       if( neset == 0 ) then
          do while( words(1) /= 'ENDEL' )
             call ecoute('mod_output_postprocess')
          end do
       else
          ivars = 0
          postp(1) % nvaes = 0
          do while( ivars < nvars )
             ivars = ivars + 1
             if( postp(1) % woese(ivars) == 'NULL ' ) then
                postp(1) % nvaes = ivars - 1
                ivars = nvars
             end if
          end do
          nvaes = 0
          do while( words(1) /= 'ENDEL' )
             call ecoute('mod_output_postprocess')
             do ivars = 1,nvars
                if( words(1) == postp(1) % woese(ivars) ) then
                   nvaes                       = nvaes + 1 
                   postp(1) % per_setse(nvaes) = ivars     
                   postp(1) % npp_setse(ivars) = 1
                   postp(1) % paese(1:5,ivars) = param(2:6)
                end if
             end do
          end do
          postp(1) % per_nvaes = nvaes
       end if

    else if( words(1) == 'BOUND' ) then
       !
       ! Boundary sets
       !
      postp(1) % npp_stepboset = getint('STEPS',1_ip,'#When to postprocess boundary sets')

       if( nbset == 0 ) then
          do while( words(1) /= 'ENDBO' )
             call ecoute('mod_output_postprocess')
          end do
       else
          ivars = 0
          postp(1) % nvabs = 0
          do while( ivars < nvars )
             ivars = ivars + 1
             if( postp(1) % wobse(ivars) == 'NULL ' ) then
                postp(1) % nvabs = ivars - 1
                ivars = nvars
             end if
          end do
          nvabs = 0                                     
          do while( words(1) /= 'ENDBO' )
             call ecoute('mod_output_postprocess')
             do ivars = 1,nvars
                if( words(1) == postp(1) % wobse(ivars) ) then
                   nvabs                       = nvabs + 1 
                   postp(1) % per_setsb(nvabs) = ivars     
                   postp(1) % npp_setsb(ivars) = 1
                   postp(1) % pabse(1:5,ivars) = param(2:6)
                end if
             end do
          end do
          postp(1) % per_nvabs = nvabs
       end if

    else if (words(1) == 'NODES' ) then
       !
       ! Node sets
       !
       postp(1) % npp_stepnoset = getint('STEPS',1_ip,'#When to postprocess node sets')

       if( nnset == 0 ) then
          do while( words(1) /= 'ENDNO' )
             call ecoute('mod_output_postprocess')
          end do
       else
          ivars = 0
          postp(1) % nvans = 0
          do while( ivars < nvars )
             ivars = ivars + 1
             if( postp(1) % wonse(ivars) == 'NULL ' ) then
                postp(1) % nvans = ivars - 1
                ivars = nvars
             end if
          end do
          do while( words(1) /= 'ENDNO' )
             call ecoute('mod_output_postprocess')
             do ivars = 1,nvars
                if( words(1) == postp(1) % wonse(ivars) ) then
                   postp(1) % npp_setsn(ivars) = 1
                   postp(1) % panse(1:5,ivars) = param(2:6)
                end if
             end do
          end do
       end if

    else if(words(1) == 'WITNE' .and. words(2) == 'GEOME' .and. .not. exists('NUMBE') ) then
       !
       ! Witness geometries
       ! 
       postp(1) % npp_stepg = 1  ! Default: postprocess at all time step
       postp(1) % npp_stepg = getint('STEPS',1_ip,'#When to postprocess witness points')
       if( nwitg == 0 ) then
          do while( words(1) /= 'ENDWI' )
             call ecoute('mod_output_postprocess')
          end do
       else
          ivarw = 0
          postp(1) % nvawg = 0
          do while( ivarw < nvarw )
             ivarw = ivarw + 1
             if( postp(1) % wowig(ivarw) == 'NULL ' ) then
                postp(1) % nvawg = ivarw - 1
                ivarw = nvarw
             end if
          end do
          do while( words(1) /= 'ENDWI' )
             call ecoute('mod_output_postprocess')
             do ivarw = 1,nvarg
                if( words(1) == postp(1) % wowig(ivarw) ) then
                   postp(1) % npp_witng(ivarw) = 1
                   if( exists('AVERA') ) then
                      postp(1) % witng_kfl_time(ivarw)  = WITNESS_AVERAGE
                      postp(1) % witng_perio_max(ivarw) = getrea('FREQU',1.0_rp,'#Frequency for averaging')
                   end if
                   if( exists('ACCUM') ) then
                      postp(1) % witng_kfl_time(ivarw)  = WITNESS_ACCUMULATE
                      postp(1) % witng_perio_max(ivarw) = getrea('FREQU',1.0_rp,'#Frequency for averaging')
                   end if
                end if
             end do
          end do
       end if

    else if (words(1) == 'WITNE' .and. .not. exists('NUMBE') ) then
       !
       ! Witness points
       !
       postp(1) % npp_stepw = 1  ! Default: postprocess at all time step
       postp(1) % npp_stepw = getint('STEPS',1_ip,'#When to postprocess witness points')
       if( nwitn == 0 ) then
          do while( words(1) /= 'ENDWI' )
             call ecoute('mod_output_postprocess')
          end do
       else
          ivarw = 0
          postp(1) % nvawi = 0
          do while( ivarw < nvarw )
             ivarw = ivarw + 1
             if( postp(1) % wowit(ivarw) == 'NULL ' ) then
                postp(1) % nvawi = ivarw - 1
                ivarw = nvarw
             end if
          end do
          do while( words(1) /= 'ENDWI' )
             call ecoute('mod_output_postprocess')
             do ivarw = 1,nvarw
                if( words(1) == postp(1) % wowit(ivarw) ) then
                   postp(1) % npp_witne(ivarw) = 1
                end if
                if( exists('AVERA') ) then
                   postp(1) % witne_kfl_time(ivarw) = WITNESS_AVERAGE
                   postp(1) % witne_period(ivarw)   = getrea('FREQU',1.0_rp,'#Frequency for averaging')
                end if
                if( exists('ACCUM') ) then
                   postp(1) % witne_kfl_time(ivarw) = WITNESS_ACCUMULATE
                   postp(1) % witne_period(ivarw)   = getrea('FREQU',1.0_rp,'#Frequency for averaging')
                end if
             end do
          end do
       end if

    else if( words(1) == 'START' ) then
       !
       ! Starting time and step of post-process
       !
       if( exists('STEP ') ) then
          postp(1) % npp_inits = getint('STEP ',0_ip,'#Initial step to start postprocess')  
          if( postp(1) % npp_inits == 0 ) postp(1) % npp_iniso = 1
       end if
       if( exists('TIME ') ) then
          postp(1) % pos_tinit = getrea('TIME ',0.0_rp,'#Initial step to start postprocess')
       end if

    end if

  end subroutine output_postprocess_read

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Allocate memory for witness
  !> @details Allocate memory for witness points
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_allocate_sets_and_witness(imodu,postp_in)

    integer(ip),           intent(in)    :: imodu
    type(typos),           intent(inout) :: postp_in

    if(maxval(postp_in % npp_setse)>0) then
       call memory_alloca(mem_modul(1:2,imodu),'VESET','mod_output_postprocess',postp_in % veset,postp_in % nvaes+1,neset)
    end if
    if(maxval(postp_in % npp_setsb)>0) then
       call memory_alloca(mem_modul(1:2,imodu),'VBSET','mod_output_postprocess',postp_in % vbset,postp_in % nvabs+1,nbset)
    end if
    if(maxval(postp_in % npp_setsn)>0) then
       call memory_alloca(mem_modul(1:2,imodu),'VNSET','mod_output_postprocess',postp_in % vnset,postp_in % nvans,nnset)
    end if
    if(maxval(postp_in % npp_witne)>0) then
       call memory_alloca(mem_modul(1:2,imodu),'WITNE','mod_output_postprocess',postp_in % witne,postp_in % nvawi,nwitn)
    end if
    if(maxval(postp_in % npp_witng)>0) then
       call memory_alloca(mem_modul(1:2,imodu),'WITNG',            'mod_output_postprocess',postp_in % witng,            postp_in % nvawg,nwitg)
       call memory_alloca(mem_modul(1:2,imodu),'WITNG_PERIO_COUNT','mod_output_postprocess',postp_in % witng_perio_count,postp_in % nvawg,nwitg)
       call memory_alloca(mem_modul(1:2,imodu),'WITNG_DENOM',      'mod_output_postprocess',postp_in % witng_denom,      postp_in % nvawg,nwitg)
       call memory_alloca(mem_modul(1:2,imodu),'WITNG_DELDENOM',   'mod_output_postprocess',postp_in % witng_deldenom,   postp_in % nvawg,nwitg)
       call memory_alloca(mem_modul(1:2,imodu),'WITNG_DELPERIO',   'mod_output_postprocess',postp_in % witng_delperio,   postp_in % nvawg,nwitg)
    end if

  end subroutine output_postprocess_allocate_sets_and_witness

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Allocate memory for witness
  !> @details Allocate memory for witness points
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_header_sets_and_witness(postp_in)

    type(typos), pointer, intent(inout) :: postp_in(:)

    if(maxval(postp_in(1) % npp_setse)>0)&
         call outset(&
         -1_ip,postp_in(1) % lun_setse,postp_in(1) % nvaes,  postp_in(1) % npp_setse,postp_in(1) % woese,postp_in(1) % veset)
    if(maxval(postp_in(1) % npp_setsb)>0)&
         call outset(&
         -2_ip,postp_in(1) % lun_setsb,postp_in(1) % nvabs,  postp_in(1) % npp_setsb,postp_in(1) % wobse,postp_in(1) % vbset)
    if(maxval(postp_in(1) % npp_setsn)>0)&
         call outset(&
         -3_ip,postp_in(1) % lun_setsn,postp_in(1) % nvans-1_ip,postp_in(1) % npp_setsn,postp_in(1) % wonse,postp_in(1) % vnset)
    if(maxval(postp_in(1) % npp_witne)>0)&
         call outset(&
         -5_ip,postp_in(1) % lun_witne,postp_in(1) % nvawi-1_ip,postp_in(1) % npp_witne,postp_in(1) % wowit,postp_in(1) % witne)
    if(maxval(postp_in(1) % npp_witng)>0)&
         call outset(&
         -6_ip,postp_in(1) % lun_witng,postp_in(1) % nvawg-1_ip,postp_in(1) % npp_witng,postp_in(1) % wowig,postp_in(1) % witng)

  end subroutine output_postprocess_header_sets_and_witness

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Output witness
  !> @details Output witness points values
  !  used by  alya/main/Output 
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_output_witness(postp_in)

    type(typos), pointer, intent(inout) :: postp_in(:)

    if( INOTSLAVE ) then
       !
       ! Element sets
       !
       if( maxval(postp_in(1) % npp_setse) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
         if(( mod(ittim, postp_in(1) % npp_stepelset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
            call outset(&
               1_ip, postp_in(1) % lun_setse,postp_in(1) % nvaes,postp_in(1) % npp_setse,postp_in(1) % woese,postp_in(1) % veset)     
            call outset(&
               10_ip,postp_in(1) % lun_setse,postp_in(1) % nvaes,postp_in(1) % npp_setse,postp_in(1) % woese,postp_in(1) % veset)
         end if
       end if
       !
       ! Boundary sets
       !
       if( maxval(postp_in(1) % npp_setsb) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
         if(( mod(ittim, postp_in(1) % npp_stepboset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
            call outset(&
               2_ip, postp_in(1) % lun_setsb,postp_in(1) % nvabs,postp_in(1) % npp_setsb,postp_in(1) % wobse,postp_in(1) % vbset)
            call outset(&
               20_ip,postp_in(1) % lun_setsb,postp_in(1) % nvabs,postp_in(1) % npp_setsb,postp_in(1) % wobse,postp_in(1) % vbset)
         end if
       end if
       !
       ! Node sets
       !
       if( maxval(postp_in(1) % npp_setsn) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
         if(( mod(ittim, postp_in(1) % npp_stepnoset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
            call outset(&
               3_ip, postp_in(1) % lun_setsn,postp_in(1) % nvans-1_ip,postp_in(1) % npp_setsn,postp_in(1) % wonse,postp_in(1) % vnset)           
            call outset(&
               30_ip,postp_in(1) % lun_setsn,postp_in(1) % nvans-1_ip,postp_in(1) % npp_setsn,postp_in(1) % wonse,postp_in(1) % vnset)
         end if
       end if
       !
       ! Witness points
       !
       if( maxval(postp_in(1) % npp_witne) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
          if(( mod(ittim, postp_in(1) % npp_stepw) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
             call outset(& ! writes header
                  5_ip, postp_in(1) % lun_witne,postp_in(1) % nvawi-1_ip,postp_in(1) % npp_witne,postp_in(1) % wowit,postp_in(1) % witne)           
             call witness_point_averaging_ini()
             call outset(& ! writes results
                  50_ip,postp_in(1) % lun_witne,postp_in(1) % nvawi-1_ip,postp_in(1) % npp_witne,postp_in(1) % wowit,postp_in(1) % witne)
             call witness_point_averaging_end()
          end if

       end if
       !
       ! Witness geometries
       !               
       if( maxval(postp_in(1) % npp_witng) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
          if(( mod(ittim, postp_in(1) % npp_stepg) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
             call outset(&
                  6_ip, postp_in(1) % lun_witng,postp_in(1) % nvawg-1_ip,postp_in(1) % npp_witng,postp_in(1) % wowig,postp_in(1) % witng)
             call witness_geometry_averaging_ini()
             call outset(&
                  60_ip,postp_in(1) % lun_witng,postp_in(1) % nvawg-1_ip,postp_in(1) % npp_witng,postp_in(1) % wowig,postp_in(1) % witng)
             call witness_geometry_averaging_end()

          end if
       end if

    end if

  end subroutine output_postprocess_output_witness

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Element set parallelization
  !> @details Element set parallelization
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_element_sets_parall()

    integer(ip)         :: nvaes,ivars 
    real(rp),   pointer :: per_veset(:,:)

    if( IPARALL .and. maxval(postp(1) % npp_setse) > 0 ) then

       nullify(per_veset)
       allocate(per_veset(postp(1) % per_nvaes+1,neset))
       do nvaes = 1,postp(1) % per_nvaes
          ivars = postp(1) % per_setse(nvaes)
          per_veset(nvaes,:) = postp(1) % veset(ivars,:) 
       end do
       per_veset(postp(1) % per_nvaes+1,:) = postp(1) % veset(postp(1) % nvaes+1,:)

       call PAR_SUM(postp(1) % per_nvaes+1_ip,neset,per_veset,'IN MY CODE') 

       do nvaes = 1,postp(1) % per_nvaes
          ivars = postp(1) % per_setse(nvaes)
          postp(1) % veset(ivars,:) = per_veset(nvaes,:) 
       end do
       postp(1) % veset(postp(1) % nvaes+1,:) = per_veset(postp(1) % per_nvaes+1,:) 
       deallocate(per_veset)

    end if

  end subroutine output_postprocess_element_sets_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Element set parallelization
  !> @details Element set parallelization
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_boundary_sets_parall()

    integer(ip)         :: nvabs,ivars 
    real(rp),   pointer :: per_vbset(:,:)

    if( IPARALL .and. maxval(postp(1) % npp_setsb) > 0 ) then

       nullify(per_vbset)
       allocate(per_vbset(postp(1) % per_nvabs+1,nbset))
       do nvabs = 1,postp(1) % per_nvabs
          ivars = postp(1) % per_setsb(nvabs)
          per_vbset(nvabs,:) = postp(1) % vbset(ivars,:) 
       end do
       per_vbset(postp(1) % per_nvabs+1,:) = postp(1) % vbset(postp(1) % nvabs+1,:)

       call PAR_SUM(postp(1) % per_nvabs+1_ip,nbset,per_vbset,'IN MY CODE')

       do nvabs = 1,postp(1) % per_nvabs
          ivars = postp(1) % per_setsb(nvabs)
          postp(1) % vbset(ivars,:) = per_vbset(nvabs,:) 
       end do
       postp(1) % vbset(postp(1) % nvabs+1,:) = per_vbset(postp(1) % per_nvabs+1,:) 
       deallocate(per_vbset)

    end if

  end subroutine output_postprocess_boundary_sets_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Element set parallelization
  !> @details Element set parallelization
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_node_sets_parall()

    if( IPARALL .and. maxval(postp(1) % npp_setsn) > 0 ) then
       nparr =  postp(1) % nvans*nnset
       vnset => postp(1) % vnset
       call par_comset(1_ip)   
    end if

  end subroutine output_postprocess_node_sets_parall


  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Witness point parallelization
  !> @details Witness point parallelization
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_witness_parall(postp_in)

    type(typos), pointer, intent(inout) :: postp_in(:)

    postp => postp_in
    witne => postp_in(1) % witne

    if( IPARALL .and. nwitn > 0 .and. maxval(postp_in(1) % npp_witne) > 0 ) then
       call par_comset(2_ip)  
    end if

  end subroutine output_postprocess_witness_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a variable is postprocessed
  !> @details Check if a variable is postprocessed at one moment
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_cancel_variable_postprocess(ivara)

    integer(ip), intent(in) :: ivara

    postp(1) % npp_stepi(ivara,:)   = 0
    postp(1) % pos_times(:,ivara,:) = 0.0_rp
    postp(1) % pos_perio(ivara,:)   = 0.0_rp

  end subroutine output_postprocess_cancel_variable_postprocess

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Sets reordering
  !> @details Redefine sets to be postprocessed... just in case modules have
  !>          redefined NPP_SETSB outside of this subroutine
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_sets_ordering()

    integer(ip) :: ivars,nvabs,nvaes

    nvaes = 0                
    do ivars = 1,nvars
       if( postp(1) % npp_setse(ivars) == 1 ) then
          nvaes                       = nvaes + 1
          postp(1) % per_setse(nvaes) = ivars    
       end if
    end do
    postp(1) % per_nvaes = nvaes

    nvabs = 0                
    do ivars = 1,nvars
       if( postp(1) % npp_setsb(ivars) == 1 ) then
          nvabs                       = nvabs + 1 
          postp(1) % per_setsb(nvabs) = ivars     
       end if
    end do
    postp(1) % per_nvabs = nvabs

  end subroutine output_postprocess_sets_ordering

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a witness point is output
  !> @details Check if a witness point is output
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function output_postprocess_check_witness_output()

    output_postprocess_check_witness_output = .false.
    if( maxval(postp(1) % npp_witne) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
       if(( mod(ittim, postp(1) % npp_stepw) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
          output_postprocess_check_witness_output = .true.
       end if
    end if

  end function output_postprocess_check_witness_output

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a variable should be postprocessed
  !> @details Check if a variable should be postprocessed
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function output_postprocess_check_variable_postprocess_now(variable_number,FORCE_POSTPROCESS,MESH_ID)

    integer(ip), intent(in)           :: variable_number
    logical(lg), intent(in), optional :: FORCE_POSTPROCESS
    integer(ip), intent(in), optional :: MESH_ID
    integer(ip)                       :: ivari,itime
    integer(ip)                       :: ivara,imesh
    integer(ip)                       :: imesh_ini,imesh_end
    logical(lg)                       :: if_force_postprocess

    if( present(FORCE_POSTPROCESS) ) then
       if_force_postprocess = FORCE_POSTPROCESS
    else
       if_force_postprocess = .false.
    end if
    if( present(MESH_ID) ) then
       imesh_ini = MESH_ID
       imesh_end = MESH_ID
    else
       imesh_ini = 0
       imesh_end = nwith
    end if

    ivara     = variable_number
    ivari     = variable_number
    ivara     = 0
    kfl_ivari = 0

    if( if_force_postprocess ) then
       ivara = ivari
    else
       do imesh = imesh_ini,imesh_end
          if( postp(1) % pos_alrea(ivari,imesh) /= 2 ) then

             if( ittyp == ITASK_INITIA .and. kfl_rstar < 2 ) then  ! outputs if INITI or no restart
                !
                ! Initial solution
                !
                if( postp(1) % npp_iniso == 1 .and. postp(1) % npp_stepi(ivari,imesh) > 0 ) then  
                   ivara = ivari
                   postp(1) % pos_alrea(ivari,imesh) = 1 ! To avoid repostprocess when no time step is required
                end if

             else if( ittyp == ITASK_ENDTIM ) then
                !
                ! End of a time step
                !
                if( cutim >= postp(1) % pos_tinit ) then
                   postp(1) % pos_alrea(ivari,imesh) = 0
                   !
                   ! At a given time step
                   !           
                   if( &
                        ittim >= postp(1) % npp_inits .and. &
                        postp(1) % pos_alrea(ivari,imesh) == 0 .and. &
                        postp(1) % npp_stepi(ivari,imesh) > 0 ) then     
                      if( mod(ittim, postp(1) % npp_stepi(ivari,imesh)) == 0 ) then
                         ivara = ivari
                         postp(1) % pos_alrea(ivari,imesh) = 1
                      end if
                   end if
                   !
                   ! At a given time period (this check should be before the posrt at a given time)
                   !                   
                   if( ittim == 1 ) postp(1) % pos_times(1,ivari,imesh) = postp(1) % pos_perio(ivari,imesh)                   
                   if( cutim >= postp(1) % pos_tinit .and. postp(1) % pos_alrea(ivari,imesh) == 0 ) then
                      if(    abs(postp(1) % pos_times(1,ivari,imesh)-cutim) < (0.6_rp*dtime).and.&
                           &     postp(1) % pos_perio(ivari,imesh)          > 0.0_rp) then
                         ivara = ivari
                         postp(1) % pos_alrea(ivari,imesh)   = 1 
                         postp(1) % pos_times(1,ivari,imesh) = postp(1) % pos_times(1,ivari,imesh) + postp(1) % pos_perio(ivari,imesh) 
                     end if
                   end if
                   !
                   ! At a given time
                   !
                   if( postp(1) % pos_alrea(ivari,imesh) == 0 ) then  
                      do itime = 1,nvart
                         if(   abs( postp(1) % pos_times(itime,ivari,imesh)-cutim) < (0.5_rp*dtime) .and. &
                              &     postp(1) % pos_times(itime,ivari,imesh)        > 0.0_rp) then
                            ivara = ivari
                            postp(1) % pos_alrea(ivari,imesh) = 1
                         end if
                      end do
                   end if

                end if

             else if( ittyp == ITASK_ENDRUN ) then
                !
                ! End of the run
                !        
                if( postp(1) % npp_stepi(ivari,imesh) /= 0 .and. postp(1) % pos_alrea(ivari,imesh) == 0 .and. &
                     ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
                   ivara = ivari   
                end if

             end if
          end if

          if( ivara == ivari .and. postp(1) % kfl_oonce(ivari) == 1 ) then
             postp(1) % pos_alrea(ivari,imesh) = 2
          end if
          
       end do
    end if

    kfl_ivari(1) = ivara

    !-------------------------------------------------------------------
    !
    ! Voxel: Initial solution, end of time step, end of the run
    !
    !-------------------------------------------------------------------

    ivara = 0

    do imesh = imesh_ini,imesh_end
       if( postp(1) % vox_alrea(ivari,imesh) /= 2 ) then
          if( ittyp == ITASK_INITIA .and. kfl_rstar == 0 ) then
             !
             ! Initial solution
             !
             if( postp(1) % npp_iniso == 1 .and. postp(1) % vox_stepi(ivari,imesh) > 0 ) then  
                ivara = ivari
                postp(1) % vox_alrea(ivari,imesh) = 1 ! To avoid repostprocess when no time step is required
             end if

          else if( ittyp == ITASK_ENDTIM ) then
             !
             ! End of a time step
             !
             if( cutim >= postp(1) % pos_tinit ) then
                postp(1) % vox_alrea(ivari,imesh) = 0
                !
                ! At a given time step
                !           
                if( &
                     ittim >= postp(1) % npp_inits .and. &
                     postp(1) % vox_alrea(ivari,imesh) == 0 .and. &
                     postp(1) % vox_stepi(ivari,imesh) > 0 ) then     
                   if( mod(ittim, postp(1) % vox_stepi(ivari,imesh)) == 0 ) then
                      ivara = ivari
                      postp(1) % vox_alrea(ivari,imesh) = 1
                   end if
                end if
                !
                ! At a given time period
                !
                if( ittim == 1 ) postp(1) % vox_times(1,ivari,imesh) = postp(1) % vox_perio(ivari,imesh)
                if( cutim >= postp(1) % pos_tinit .and. postp(1) % vox_alrea(ivari,imesh) == 0 ) then  
                   if(    abs(postp(1) % vox_times(1,ivari,imesh)-cutim) < (0.6_rp*dtime).and.&
                        &     postp(1) % vox_perio(ivari,imesh)          > 0.0_rp) then
                      ivara = ivari
                      postp(1) % vox_alrea(ivari,imesh)   = 1 
                      postp(1) % vox_times(1,ivari,imesh) = postp(1) % vox_times(1,ivari,imesh) + postp(1) % vox_perio(ivari,imesh) 
                   end if
                end if
                !
                ! At a given time
                !
                if( postp(1) % vox_alrea(ivari,imesh) == 0 ) then  
                   do itime = 1,nvart
                      if(   abs( postp(1) % vox_times(itime,ivari,imesh)-cutim) < (0.5_rp*dtime) .and. &
                           &     postp(1) % vox_times(itime,ivari,imesh)        > 0.0_rp) then
                         ivara = ivari
                         postp(1) % vox_alrea(ivari,imesh) = 1
                      end if
                   end do
                end if
 
             end if

          else if( ittyp == ITASK_ENDRUN ) then
             !
             ! End of the run
             !        
             if( postp(1) % vox_stepi(ivari,imesh) /= 0 .and. postp(1) % vox_alrea(ivari,imesh) == 0 .and. &
                  ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
                ivara = ivari   
             end if

          end if
       end if

       if( ivara == ivari .and. postp(1) % kfl_oonce(ivari) == 1 ) then
          postp(1) % vox_alrea(ivari,imesh) = 2
       end if
       
    end do

    kfl_ivari(2) = ivara

    ivara = maxval(kfl_ivari)

    if( ivara == 0 ) then
       output_postprocess_check_variable_postprocess_now = .false.
    else
       output_postprocess_check_variable_postprocess_now = .true.
    end if

  end function output_postprocess_check_variable_postprocess_now

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a variable is postprocessed
  !> @details Check if a variable is postprocessed at one moment
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function output_postprocess_check_variable_postprocess(ivara,VARIABLE_NAME)

    integer(ip),      intent(in), optional :: ivara
    character(len=*), intent(in), optional :: VARIABLE_NAME
    integer(ip)                            :: jvara,imesh

    output_postprocess_check_variable_postprocess = .true.
    
    if( present(ivara) ) then
       jvara = ivara
    else if( present(VARIABLE_NAME) ) then
       jvara = 1
       do while( trim(postp(1) % wopos(1,jvara)) /= trim(VARIABLE_NAME) )
          jvara = jvara + 1
          if( jvara > size(postp(1) % wopos(1,:)) ) then
             call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME))
          end if
       end do
    else
       call runend('MOD_OUTPUT_POTPROCESS: MISSING ARGUMENT')
    end if

    do imesh = 0,nwith
       if(  postp(1) % npp_stepi(jvara,imesh)           /= 0      .or. &
            maxval(postp(1) % pos_times(:,jvara,imesh)) >  0.0_rp .or. &
            postp(1) % pos_perio(jvara,imesh)           /= 0.0_rp ) then
          return
       else
          output_postprocess_check_variable_postprocess = .false.
       end if
    end do
    
  end function output_postprocess_check_variable_postprocess
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @date    2022-02-24
  !> @brief   Check if a witness variable is postprocessed
  !> @details Check if a witness variable is postprocessed
  !> 
  !-----------------------------------------------------------------------
  
  logical(lg) function output_postprocess_check_variable_witness(ivara,VARIABLE_NAME,NCHAR)

    integer(ip),      intent(in), optional :: ivara
    character(len=*), intent(in), optional :: VARIABLE_NAME
    integer(ip),      intent(in), optional :: NCHAR
    integer(ip)                            :: jvara

    output_postprocess_check_variable_witness = .false.
    jvara = 1
    
    if( present(ivara) ) then
       jvara = ivara
    else if( present(VARIABLE_NAME) ) then
       if( present(NCHAR) ) then
          do while( trim(postp(1) % wowit(jvara)(1:NCHAR)) /= trim(VARIABLE_NAME(1:NCHAR)) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_witne) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME(1:NCHAR)))
             end if
          end do
       else
          do while( trim(postp(1) % wowit(jvara)) /= trim(VARIABLE_NAME) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_witne) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME))
             end if
          end do
       end if
    else
       call runend('MOD_OUTPUT_POTPROCESS: MISSING ARGUMENT')
    end if

    if( maxval(postp(1) % npp_witne) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
       if(( mod(ittim, postp(1) % npp_stepw) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
          if( postp(1) % npp_witne(jvara) /= 0 ) then
             output_postprocess_check_variable_witness = .true.
          end if
       end if
    end if
    
  end function output_postprocess_check_variable_witness
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @date    2022-02-24
  !> @brief   Check if a node set variable is postprocessed
  !> @details Check if a node set variable is postprocessed
  !> 
  !-----------------------------------------------------------------------
  
  logical(lg) function output_postprocess_check_variable_node_sets(ivara,VARIABLE_NAME,NCHAR)

    integer(ip),      intent(in), optional :: ivara
    character(len=*), intent(in), optional :: VARIABLE_NAME
    integer(ip),      intent(in), optional :: NCHAR
    integer(ip)                            :: jvara

    output_postprocess_check_variable_node_sets = .false.
    jvara = 1

    if( present(ivara) ) then
       jvara = ivara
    else if( present(VARIABLE_NAME) ) then
       if( present(NCHAR) ) then
          do while( trim(postp(1) % wonse(jvara)(1:NCHAR)) /= trim(VARIABLE_NAME(1:NCHAR)) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_setsn) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME(1:NCHAR)))
             end if
          end do
       else
          do while( trim(postp(1) % wonse(jvara)) /= trim(VARIABLE_NAME) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_setsn) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME))
             end if
          end do
       end if
    else
       call runend('MOD_OUTPUT_POTPROCESS: MISSING ARGUMENT')
    end if

    if( maxval(postp(1) % npp_setsn) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
       if(( mod(ittim, postp(1) % npp_stepnoset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
          if( postp(1) % npp_setsn(jvara) /= 0 ) then
             output_postprocess_check_variable_node_sets = .true.
          end if
       end if
    end if
    
  end function output_postprocess_check_variable_node_sets

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @date    2022-02-24
  !> @brief   Check if a boundary set variable is postprocessed
  !> @details Check if a boundary set variable is postprocessed
  !> 
  !-----------------------------------------------------------------------
  
  logical(lg) function output_postprocess_check_variable_boundary_sets(ivara,VARIABLE_NAME,NCHAR)

    integer(ip),      intent(in), optional :: ivara
    character(len=*), intent(in), optional :: VARIABLE_NAME
    integer(ip),      intent(in), optional :: NCHAR
    integer(ip)                            :: jvara

    output_postprocess_check_variable_boundary_sets = .false.
    jvara = 1
    
    if( present(ivara) ) then
       jvara = ivara
    else if( present(VARIABLE_NAME) ) then
       if( present(NCHAR) ) then
          do while( trim(postp(1) % wobse(jvara)(1:NCHAR)) /= trim(VARIABLE_NAME(1:NCHAR)) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_setsb) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME(1:NCHAR)))
             end if
          end do
       else
          do while( trim(postp(1) % wobse(jvara)) /= trim(VARIABLE_NAME) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_setsb) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME))
             end if
          end do
       end if
    else
       call runend('MOD_OUTPUT_POTPROCESS: MISSING ARGUMENT')
    end if
    
    if( maxval(postp(1) % npp_setsb) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
       if(( mod(ittim, postp(1) % npp_stepboset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
          if( postp(1) % npp_setsb(jvara) /= 0 ) then
             output_postprocess_check_variable_boundary_sets = .true.
          end if
       end if
    end if
    
  end function output_postprocess_check_variable_boundary_sets

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @date    2022-02-24
  !> @brief   Check if a element set variable is postprocessed
  !> @details Check if a element set variable is postprocessed
  !> 
  !-----------------------------------------------------------------------
  
  logical(lg) function output_postprocess_check_variable_element_sets(ivara,VARIABLE_NAME,NCHAR)

    integer(ip),      intent(in), optional :: ivara
    character(len=*), intent(in), optional :: VARIABLE_NAME
    integer(ip),      intent(in), optional :: NCHAR
    integer(ip)                            :: jvara

    output_postprocess_check_variable_element_sets = .false.
    jvara = 1
    
    if( present(ivara) ) then
       jvara = ivara
    else if( present(VARIABLE_NAME) ) then
       if( present(NCHAR) ) then
          do while( trim(postp(1) % woese(jvara)(1:NCHAR)) /= trim(VARIABLE_NAME(1:NCHAR)) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_setse) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME(1:NCHAR)))
             end if
          end do
       else
          do while( trim(postp(1) % woese(jvara)) /= trim(VARIABLE_NAME) )
             jvara = jvara + 1
             if( jvara > size(postp(1) % npp_setse) ) then
                call runend('MOD_OUTPUT_POTPROCESS: VARIABLE NOT FOUND '//trim(VARIABLE_NAME))
             end if
          end do
       end if
    else
       call runend('MOD_OUTPUT_POTPROCESS: MISSING ARGUMENT')
    end if

    if( maxval(postp(1) % npp_setse) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
       if(( mod(ittim, postp(1) % npp_stepelset) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
          if( postp(1) % npp_setse(jvara) /= 0 ) then
             output_postprocess_check_variable_element_sets = .true.
          end if
       end if
    end if

  end function output_postprocess_check_variable_element_sets
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Read a geometry
  !> @details Read a geometry
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_read_geometry(kfl_geometry,param_geo)

    integer(ip), intent(out) :: kfl_geometry
    real(rp),    intent(out) :: param_geo(:)
    real(rp)                 :: normal(3),rr

    if( words(1) == 'SPHER' .or. words(1) == 'CIRCL' ) then
       !
       ! Sphere or circle
       !
       kfl_geometry = WITNESS_SPHERE
       param_geo(1) = getrea('RADIU',0.0_rp,'#radius of sphere')
       param_geo(2) = getrea('X    ',0.0_rp,'#center of sphere x')
       param_geo(3) = getrea('Y    ',0.0_rp,'#center of sphere y')
       param_geo(4) = getrea('Z    ',0.0_rp,'#center of sphere z')
       
    else if( words(1) == '3DCIR' .or. words(1) == '3DRIN' ) then
       !
       ! 3D circle
       !
       if(      words(1) == '3DCIR' ) then
          kfl_geometry = WITNESS_3D_CIRCLE
       else if( words(1) == '3DRIN' ) then
          kfl_geometry = WITNESS_3D_RING          
       end if
       param_geo(1) = getrea('RADIU', 0.0_rp,'#radius of circle')
       param_geo(2) = getrea('X    ', 0.0_rp,'#center of circle x')
       param_geo(3) = getrea('Y    ', 0.0_rp,'#center of circle y')
       param_geo(4) = getrea('Z    ', 0.0_rp,'#center of circle z')
       param_geo(5) = getrea('NX   ', 0.0_rp,'#normal of circle x')
       param_geo(6) = getrea('NY   ', 0.0_rp,'#normal of circle y')
       param_geo(7) = getrea('NZ   ', 0.0_rp,'#normal of circle z')
       param_geo(8) = getrea('NR   ',10.0_rp,'#nodes in r direction')
       param_geo(9) = getrea('NTHET',10.0_rp,'#node in theta direction')
       
       call maths_normalize_vector(3_ip,param_geo(5:7),rr)
       
    else if( words(1) == 'PLANE' ) then
       !
       ! Plane
       !
       kfl_geometry = WITNESS_PLANE
       param_geo(1) = getrea('A    ',0.0_rp,'#plane equation a')
       param_geo(2) = getrea('B    ',0.0_rp,'#plane equation b')
       param_geo(3) = getrea('C    ',0.0_rp,'#plane equation c')
       param_geo(4) = getrea('D    ',0.0_rp,'#plane equation d')
       if( exists('X    ') ) then
          normal    = 0.0_rp
          normal(1) = getrea('X    ',0.0_rp,'#plane equation x')
          normal(2) = getrea('Y    ',0.0_rp,'#plane equation y')
          normal(3) = getrea('Z    ',0.0_rp,'#plane equation z')
          param_geo(4) = -normal(1)*param_geo(1)-normal(2)*param_geo(2)-normal(3)*param_geo(3)
       end if
       
    else if( words(1) == 'BAR3D' ) then
       !
       ! Bar3d
       !
       kfl_geometry = WITNESS_BAR3D
       param_geo(1) = getrea('RADIA',5.0_rp,'#Number elements in radial direction')
       param_geo(2) = getrea('CIRCU',6.0_rp,'#Number elements in circumference direction')
       
    else if( words(1) == 'BOUND' ) then
       !
       ! Boundary, coarse boundary, boundary set
       !
       kfl_geometry = WITNESS_BOUNDARY
       param_geo(1) = 0.0_rp
       if(      exists('NODES') ) then
          param_geo(1) = getrea('NODES',0.0_rp,'#Total number of nodes')          
       else if( exists('SET  ')   ) then
          kfl_geometry = WITNESS_BOUNDARY_SET
          param_geo(1) = getrea('SET  ',1.0_rp,'#set number')          
       end if
       
    else if( words(1) == 'MATER' ) then
       !
       ! Material
       !
       kfl_geometry = WITNESS_MATERIAL
       param_geo(1) = getrea('MATER',1.0_rp,'#Material set')
       
    else if( words(1) == 'ELEME' ) then
       !
       ! Element set
       !
       kfl_geometry = WITNESS_ELEMENT_SET
       if( exists('SET  ') ) then
          param_geo(1) = getrea('SET  ',1.0_rp,'#Element set')
       else
          param_geo(1) = getrea('ELEME',1.0_rp,'#Element set')
       end if
       
    else if( words(1) == 'BOX  ' ) then
       !
       ! Box
       !
       kfl_geometry = WITNESS_BOX
       param_geo(1) = getrea('XMIN ',0.0_rp,'#lower limit in x')
       param_geo(1+ndime)= getrea('XMAX ',0.0_rp,'#upper limit in x')
       if ( ndime >=2 ) then
          param_geo(2)       = getrea('YMIN ',0.0_rp,'#lower limit in y')
          param_geo(2+ndime) = getrea('YMAX ',0.0_rp,'#upper limit in y')
       end if
       if ( ndime >=3 ) then
          param_geo(3)       = getrea('ZMIN ',0.0_rp,'#lower limit in z')
          param_geo(3+ndime) = getrea('ZMAX ',0.0_rp,'#upper limit in z')
       end if
       param_geo(7)       = getrea('NX   ',0.0_rp,'#Number of boxes in x')
       param_geo(8)       = getrea('NY   ',0.0_rp,'#Number of boxes in y')
       param_geo(9)       = getrea('NZ   ',0.0_rp,'#Number of boxes in z')
       
    else if( words(1) == 'INBOX' ) then
       !
       ! In a Box
       !
       kfl_geometry = WITNESS_IN_BOX
       param_geo(1) = getrea('XMIN ',0.0_rp,'#lower limit in x')
       param_geo(1+ndime)= getrea('XMAX ',0.0_rp,'#upper limit in x')
       if ( ndime >=2 ) then
          param_geo(2)       = getrea('YMIN ',0.0_rp,'#lower limit in y')
          param_geo(2+ndime) = getrea('YMAX ',0.0_rp,'#upper limit in y')
       end if
       if ( ndime >=3 ) then
          param_geo(3)       = getrea('ZMIN ',0.0_rp,'#lower limit in z')
          param_geo(3+ndime) = getrea('ZMAX ',0.0_rp,'#upper limit in z')
       end if
       
    else if( words(1) == 'RING ' ) then
       !
       ! Ring
       !
       kfl_geometry = WITNESS_RING
       param_geo(1) = getrea('X    ',0.0_rp,'#center of ring x')
       param_geo(2) = getrea('Y    ',0.0_rp,'#center of ring y')
       param_geo(3) = getrea('Z    ',0.0_rp,'#center of ring z')
       param_geo(4) = getrea('NX   ',1.0_rp,'#normal of ring x')
       param_geo(5) = getrea('NY   ',0.0_rp,'#normal of ring y')
       param_geo(6) = getrea('NZ   ',0.0_rp,'#normal of ring z')
       param_geo(7) = getrea('R    ',0.0_rp,'#radius of ring')
       param_geo(8) = getrea('DR   ',0.0_rp,'#thickness in radial direction')
       param_geo(9) = getrea('DN   ',0.0_rp,'#thickness in normal direction')
       !
       ! Normalize normal
       !
       normal(1:ndime)      = param_geo(4:3+ndime)
       rr                   = sqrt(dot_product(normal,normal))
       param_geo(4:3+ndime) = normal(1:ndime)/rr
       
    end if

  end subroutine output_postprocess_read_geometry

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-08-27
  !> @brief   Output variables on meshes
  !> @details Output variables on meshes
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_variables(subru)

    procedure(xxx_outvar) :: subru
    integer(ip)           :: ivari,imesh

    do ivari = 1,nvarp
       do imesh = 0,nwith
          if( output_postprocess_check_variable_postprocess_now(ivari,MESH_ID=imesh) ) then
             call subru(ivari,imesh)
          end if
       end do
    end do

  end subroutine output_postprocess_variables
  
end module mod_output_postprocess
!> @}
