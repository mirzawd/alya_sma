!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!--------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_tem_memory.f90
!> @author  Adria Quintanas-Corominas
!> @author  David Oks 
!> @author  Guillaume Houzeaux
!> @date    June, 2020
!> @brief   
!> @details 
!--------------------------------------------------------------------

module mod_tubes
  !------------------------
  use def_kintyp_basic,          only : ip, rp, lg, i2p
  use def_master,                only : mmodu
  use def_master,                only : modul
  use def_master,                only : mem_modul
  use def_master,                only : kfl_modul
  use def_master,                only : INOTMASTER,kfl_paral
  use def_master,                only : ISLAVE
  use def_domain,                only : mcodb
  use mod_communications_global, only : PAR_BROADCAST
  use mod_communications_global, only : PAR_SUM
  use mod_communications_global, only : PAR_MAX
  use mod_communications,        only : PAR_SEND_RECEIVE_TO_ALL
  use mod_memory_basic,          only : memory_alloca
  use mod_ecoute,                only : ecoute
  use mod_messages,              only : messages_live
  use mod_strings,               only : integer_to_string
  use mod_strings,               only : real_to_string
  use def_inpout
  !------------------------
  implicit none
  private

  !------------------------
  
  integer(ip),         parameter :: len_file = 200
  !------------------------
  type :: typ_list_values
     real(rp), pointer           :: tubes_in(:)
     real(rp), pointer           :: tubes_out(:)
  end type typ_list_values

  type :: typ_tubes
     integer(ip)                 :: mpi
     character(len_file)         :: casename 
     character(len_file)         :: output 
   contains
     procedure,             pass :: init
     procedure,             pass :: read_data
     procedure,             pass :: parall
  end type typ_tubes
  !------------------------
  integer(ip)                    :: num_tubes
  integer(ip)                    :: kfl_exists_tubes
  integer(ip)                    :: ndifferential_timesteps
  integer(ip)                    :: num_codes(0:mmodu)
  type(i2p)                      :: list_codes(0:mmodu)
  type(typ_list_values)          :: list_values(0:mmodu)
  type(typ_tubes),       pointer :: tubes(:)
  !------------------------
  public                         :: unnitt_tubes
  public                         :: tubes_initialization
  public                         :: tubes_register
  public                         :: tubes_parall
  public                         :: tubes_recv_value_from_module
  public                         :: tubes_give_value_to_module
  public                         :: tubes_total_number
  public                         :: tubes_number
  public                         :: tubes_doiter
  public                         :: tubes_iniunk
  public                         :: tubes_readat
  !------------------------
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Initialize a tube
  !> @details Initialization of a tube
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init(tube)

    class(typ_tubes) :: tube

    tube % mpi      = 0
    tube % casename = ''
    
  end subroutine init
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Broadcast tube info
  !> @details Broadcast tube info to slaves
  !> 
  !-----------------------------------------------------------------------
  
  subroutine parall(tube)

    class(typ_tubes) :: tube

    call PAR_BROADCAST(tube % mpi)
    call PAR_BROADCAST(len_file,tube % casename,'IN MY CODE')
    
  end subroutine parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Read tube data
  !> @details Read tube data
  !> 
  !-----------------------------------------------------------------------
  
  subroutine read_data(tube)

    class(typ_tubes) :: tube

    if(      words(1) == 'MPI  ' ) then
       tube % mpi = getint('MPI  ',1_ip,'#MPI number')
    else if( words(1) == 'CASEN' ) then
       tube % casename = trim(getcha('CASEN','NULL ','#Case name'))
    end if
    
  end subroutine read_data
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Read data
  !> @details Read tubes configuration data
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tubes_readat()

    integer(ip) :: itube

    if( words(2) /= 'OFF  ' ) then

       if( exists('NUMBE') ) then
          num_tubes = getint('NUMBE',0_ip,'#Number of tubes')
       else
          num_tubes = getint('TUBES',0_ip,'#Number of tubes')
       end if
       if( num_tubes > 0 ) then
          kfl_exists_tubes = 1
          allocate( tubes(num_tubes) )
          do itube = 1,num_tubes
             call tubes(itube) % init()
          end do

          do while( words(1) /= 'ENDTU' )

             if( words(1) == 'TUBE ' ) then
                if( exists('NUMBE') ) then
                   itube = getint('NUMBE',1_ip,'#Tube number')
                else
                   itube = getint('TUBE ',1_ip,'#Tube number')
                end if
                if( itube == 0 ) call runend('MOD_TUBES: TUBE NUMBER NOT SPECIFIED')
                do while( words(1) /= 'ENDTU' )
                   call tubes(itube) % read_data()
                   call ecoute('mod_tubes')
                end do
             end if

             call ecoute('mod_tubes')
          end do
       end if

    end if

  end subroutine tubes_readat
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   doiter
  !> @details Solve tubes
  !> 
  !-----------------------------------------------------------------------

  subroutine tubes_doiter()
    
    integer(ip)     :: imodu,itime,ii,itube,isolv
!    integer(ip)     :: num_total_tubes
    real(rp)        :: dQ, Qiter,P
    real(rp), save  :: Qold=0.0_rp

    if( kfl_exists_tubes /= 0 .and. num_tubes > 0 ) then
       call messages_live('SOLVE TUBES','START SECTION')
       do imodu = 0,mmodu
          if( kfl_modul(imodu) /= 0 .and. num_codes(imodu) > 0 ) then
             do ii = 1,num_codes(imodu)   !num_total_tubes
                itube = list_codes(imodu) % l(1,ii)
                isolv = 0
                if( tubes(itube) % mpi == kfl_paral ) then
                   isolv = kfl_paral
                else
                   isolv = 0
                end if
                call PAR_SUM(isolv)
                call messages_live('SOLVE TUBE '//integer_to_string(itube)//' BY MPI RANK '//integer_to_string(itube))
                if( tubes(itube) % mpi == kfl_paral ) then
                   ! Runs local timeloop for diferential timestepping 
                   dQ    = (list_values(imodu) % tubes_in(ii) - Qold) / real(ndifferential_timesteps,rp)
                   Qiter = Qold
                   do itime = 1,ndifferential_timesteps
                      Qiter = Qiter + dQ
                      P     = 0.0_rp
#ifdef TUBES
                      call advance_timestep_tubes(Qiter,P)
#endif
                      list_values(imodu) % tubes_out(ii) = P
                    end do
                   Qold = Qiter
                end if
             end do
             ! solve tubes for module imodu
             ! Output es list_values(imodu) % tubes_out
             ! Your tube number is list_codes(imodu) % l(1,itube)
          end if
       end do
       call messages_live('SOLVE TUBES','END SECTION')
    end if
    
  end subroutine tubes_doiter
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Tube number
  !> @details Give a tube number according to a boundary code
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function tubes_number(icode_module)

    integer(ip), intent(in)  :: icode_module
    integer(ip)              :: ii(1)
    integer(ip)              :: imodu
    
    imodu = modul
    tubes_number = 0
    if( num_codes(imodu) <= 0 ) then
       call runend('MOD_TUBES: TUBES WAS NOT REGISTERED FOR THIS MODULE 1')
    else
       if (.not.(any(list_codes(imodu) % l(2,:)==icode_module))) return
       ii = minloc(list_codes(imodu) % l(2,:),list_codes(imodu) % l(2,:)==icode_module)
       if( ii(1) <= 0 ) then
          return
       else
          tubes_number = ii(1) 
       end if
    end if

  end function tubes_number
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Total numbner of tubes
  !> @details Total numbner of tubes
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function tubes_total_number()
    tubes_total_number = num_codes(modul)
    
  end function tubes_total_number
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Give value from module
  !> @details Give a value from a module
  !> 
  !-----------------------------------------------------------------------
  
  function tubes_give_value_to_module(icode_module) result(xvalu)

    integer(ip), intent(in)  :: icode_module
    real(rp)                 :: xvalu
    integer(ip)              :: ii(1)
    integer(ip)              :: imodu

    imodu = modul
    if( num_codes(imodu) <= 0 ) then
       call runend('MOD_TUBES: TUBES WAS NOT REGISTERED FOR THIS MODULE 1')
    else
       ii = minloc(list_codes(imodu) % l(2,:),list_codes(imodu) % l(2,:)==icode_module)
       if( ii(1) <= 0 ) then
          call runend('MOD_TUBES: THIS CODE WAS NOT FOUND')
       else
          xvalu = list_values(imodu) % tubes_out(ii(1))
       end if
    end if
    
  end function tubes_give_value_to_module
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Receive value from module
  !> @details Receive a value from a module
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tubes_recv_value_from_module(xvalu)

    real(rp),    intent(in) :: xvalu(:)
    integer(ip)             :: imodu,itube

    imodu = modul

    do itube = 1,size(xvalu)
       list_values(imodu) % tubes_in(itube) = xvalu(itube)
    end do
    
  end subroutine tubes_recv_value_from_module
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Parall
  !> @details Master send results to slaves
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tubes_parall()

    integer(ip) :: imodu,itube

    call PAR_BROADCAST(kfl_exists_tubes)

    if( kfl_exists_tubes /= 0 ) then
       call PAR_BROADCAST(mmodu+1_ip,num_codes)
       call PAR_BROADCAST(num_tubes)       
       do imodu = 0,mmodu
          if( num_codes(imodu) > 0 ) then
             if(  .not. associated(list_codes(imodu)%l) ) then
                call memory_alloca(mem_modul(1:2,modul),'LIST_CODES % L','tubes_parall',list_codes(imodu) % l,2_ip,num_codes(imodu))
             end if
             call PAR_BROADCAST(2_ip,num_codes(imodu),list_codes(imodu)%l)
          end if
       end do
       
       if( ISLAVE ) allocate( tubes(num_tubes) )
       do itube = 1,num_tubes
          call tubes(itube) % parall()
       end do
    end if
    
  end subroutine tubes_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Tubes allocation 
  !> @details Allocate memory for tubes in/out values
  !> 
  !-----------------------------------------------------------------------

  subroutine tubes_allocate_inout_values()

    integer(ip) :: imodu
    
    do imodu = 0,mmodu
       if( num_codes(imodu) > 0 ) then
          call memory_alloca(mem_modul(1:2,modul),'LIST_VALUES % TUBES_IN' ,'tubes_parall',list_values(imodu) % tubes_in ,num_codes(imodu))
          call memory_alloca(mem_modul(1:2,modul),'LIST_VALUES % TUBES_OUT','tubes_parall',list_values(imodu) % tubes_out,num_codes(imodu))
       end if
    end do

  end subroutine tubes_allocate_inout_values
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Register
  !> @details Register codes 
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tubes_register(icode_tubes,icode_module,imodu)

    integer(ip), intent(in) :: icode_tubes
    integer(ip), intent(in) :: icode_module
    integer(ip), intent(in) :: imodu

    if( .not. associated(list_codes(imodu)%l) ) then
       call memory_alloca(mem_modul(1:2,modul),'LIST_CODES % L','tubes_parall',list_codes(imodu) % l,2_ip,mcodb)      
    end if
    num_codes(imodu)                          = num_codes(imodu) + 1    
    kfl_exists_tubes                          = 1
    list_codes(imodu) % l(1,num_codes(imodu)) = icode_tubes
    list_codes(imodu) % l(2,num_codes(imodu)) = icode_module
    
  end subroutine tubes_register
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-16
  !> @brief   Initialization
  !> @details Initialization of tubes variables
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tubes_initialization()

    integer(ip) :: imodu
    
    do imodu = 0,mmodu
       nullify(list_codes(imodu)%l)
       nullify(list_values(imodu)%tubes_in)
       nullify(list_values(imodu)%tubes_out)
    end do
    num_codes        = 0
    kfl_exists_tubes = 0
    nullify(tubes)
    
  end subroutine tubes_initialization
  
  !--------------------------------------------------------------------
  !> @author  David Oks 
  !> @author  Adria Quintanas-Corominas
  !> @date    June, 2020
  !> @brief   Initialises Tubes and computes it's timestep.
  !> @details 
  !--------------------------------------------------------------------
  
  subroutine tubes_iniunk()

    !------------------------
    implicit none
    !------------------------
#ifdef TUBES
    integer(ip)         :: imodu, itube
    real(rp)            :: dt_tubes
    real(rp),   pointer :: dt_all(:,:)
    !------------------------

    if( num_tubes > 0 ) then

       nullify(dt_all)
       allocate(dt_all(2,num_tubes))
       dt_all = 0.0_rp
       call messages_live('TUBES','START SECTION')
       !
       ! Allocate memory
       ! 
       call tubes_allocate_inout_values()

       do itube = 1,num_tubes
          call messages_live('READING FILES OF TUBE '//integer_to_string(itube)//' BY MPI RANK '//integer_to_string(tubes(itube) % mpi))        
          if( tubes(itube) % mpi == kfl_paral ) then
             call tubes_files(tubes(itube) % casename,54_4,55_4,56_4,57_4,58_4,59_4,60_4)
             ! Read original timestep from Tubes' inputs
             call get_timestep_tubes(dt_tubes)
             ! Compute timestep for Tubes
             call compute_timestep_for_tubes(itube,dt_tubes,ndifferential_timesteps,dt_all)
             ! Initialize Tubes program
             call initialize_tubes(dt_tubes)
          end if
       end do
       call PAR_MAX(dt_all)
       do itube = 1,num_tubes
          call messages_live('TUBE '//integer_to_string(itube),'START SECTION')
          call messages_live('INITIAL    TIME STEP= '//real_to_string(dt_all(1,itube),REAL_FORMAT='(e11.4)'))
          call messages_live('CALCULATED TIME STEP= '//real_to_string(dt_all(2,itube),REAL_FORMAT='(e11.4)'))          
          call messages_live('TUBE '//integer_to_string(itube),'END SECTION')
       end do
       call messages_live('END TUBES','END SECTION')
       deallocate(dt_all)
    end if


#endif

  end subroutine tubes_iniunk

  !--------------------------------------------------------------------
  !> @author  David Oks 
  !> @author  Adria Quintanas-Corominas
  !> @date    June, 2020
  !> @brief   Computes timestep for Tubes used for differential timestepping.
  !> @details 
  !--------------------------------------------------------------------
  subroutine compute_timestep_for_tubes(itube,dt_tubes, ndifferential_timesteps,dt_all)
    use def_master,       only : dtime
    !------------------------
    implicit none
    !------------------------
    integer(ip), intent(in)    :: itube
    real(rp),    intent(inout) :: dt_tubes
    integer(ip), intent(out)   :: ndifferential_timesteps
    real(rp),    intent(inout) :: dt_all(:,:)
    !------------------------

    dt_all(1,itube) = dt_tubes

    ! Calculate factor of dt_alya nearest to original dt_tubes (round downards)
    if (dtime > dt_tubes) then
       ndifferential_timesteps = int(ceiling(dtime / dt_tubes))
       dt_tubes                = dtime / real(ndifferential_timesteps, rp)
    else
       dt_tubes = dtime
    end if
    
    dt_all(2,itube) = dt_tubes
    
  end subroutine compute_timestep_for_tubes


  !--------------------------------------------------------------
  !> @author  Adria Quintanas-Corominas
  !> @author  David Oks 
  !> @date    June, 2020
  !> @brief  
  !> @details 
  !--------------------------------------------------------------
  subroutine unnitt_tubes( &
       icase  )
    !------------------------
    implicit none
    !------------------------
    integer(ip), intent(in) :: icase
#ifdef TUBES
    integer(ip)             :: i, ntmax
    !------------------------
    real(rp)                :: Q, P, t, dt, f
#endif
    !------------------------
    integer(ip), parameter  :: SINUSOIDE_INTERPOLATION = 0_ip
    real(rp),    parameter  :: pi = 4.0_rp*atan(1.0_rp) 
    !------------------------
#ifdef TUBES
    select case(icase)

    case(SINUSOIDE_INTERPOLATION)

       ! Parameters 
       f     = 5.0_rp    ! frequency [hz]
       dt    = 5e-6_rp   ! timestep [s]
       ntmax = 10000_ip ! max nmbr of timesteps

       ! run tubes externally, set a sinusoidal flowrate curve
       Q = 0.0_rp
       t = 0.0_rp
       call initialize_tubes()
       do i = 1, ntmax
          call advance_timestep_tubes(Q, P)
          t = t + dt
          Q = sin(2.0_rp*pi*f*t)
       end do

    end select
#else
    call runend('mod_tubes: Tubes not defined!')
#endif
    !------------------------
  end subroutine unnitt_tubes

end module mod_tubes
