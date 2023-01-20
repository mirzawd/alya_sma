!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @author  J.C. Cajas
!> @date    17/02/2014
!> @brief   Split MPI_COMM_WORLD in code communicators
!> @details Split MPI_COMM_WORLD in code communicators, based on the 
!>          name of the executable the code communicator is stored 
!>          to PAR_COMM_MY_CODE_ARRAY(:). This is done before knowing 
!>          the number of zones and subdomains.
!>
!>                        MPI_COMM_WORLD = PAR_COMM_UNIVERSE              |
!>                                 (THE UNIVERSE)                         |
!>                                       ||                               |
!>                                       \/                               |
!>                                                                        |
!>                                 MPI_COMM_SPLIT                         |
!>                                                                        | 
!>                                //            \\                        |        PERFORMED IN 
!>                               //              \\                       |  par_code_split_universe
!>                                                                        | 
!>                   PAR_COMM_WORLD             REST OF THE WORLD         |
!>              (ALYA's COMMUNICATOR)        (OTHER CODES COMMUNICATOR)   |
!>                        ||                                              | 
!>                        \/                                              |
!> ----------------------------------------------------------------------------------------
!>                  MPI_COMM_SPLIT                                        |
!>                                                                        |
!>                        ||                                              |
!>                        \/                                              |
!>                                                                        |        PERFORMED IN 
!> PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD = PAR_COMM_MY_CODE          |   par_code_split_world
!>                 (par_sencom.f90)                                       |
!>                                                                        |
!>                        ||  MPI_COMM_SPLIT                              |
!>                        \/                                              |
!>                                                                        |
!> PAR_COMM_COLOR_ARRAY(:) % PAR_COMM_WORLD (par_color_communicators.f90) |
!>                                                                        |
!>                                                                        |       
!> @} 
!
!----------------------------------------------------------------------

subroutine par_code_split_world()
   use def_kintyp,         only : ip
   use def_master,         only : IPARALL, intost
   use mod_parall,         only : PAR_COMM_UNIVERSE
   use mod_parall,         only : PAR_UNIVERSE_SIZE
   use mod_parall,         only : PAR_MY_UNIVERSE_RANK
   use mod_parall,         only : PAR_COMM_MY_CODE
   use mod_parall,         only : PAR_COMM_MY_CODE_WM
   use mod_parall,         only : PAR_COMM_SHARED_WM
   use mod_parall,         only : PAR_COMM_WORLD
   use mod_parall,         only : PAR_WORLD_SIZE
   use mod_parall,         only : PAR_MY_WORLD_RANK
   use mod_parall,         only : PAR_MY_CODE_RANK_WM
   use mod_parall,         only : PAR_SHARED_RANK_WM
   use mod_parall,         only : PAR_COMM_CURRENT
   use mod_parall,         only : mcode
   use mod_communications, only : PAR_INIT
   use mod_communications, only : PAR_COMM_RANK_AND_SIZE
   use mod_communications, only : PAR_GATHER
   use mod_communications, only : PAR_COMM_SPLIT
   use mod_communications, only : PAR_SCATTER
   use mod_communications, only : PAR_BROADCAST
   use mod_communications, only : PAR_DEFINE_COMMUNICATOR
   use mod_strings,        only : string_sort
   use mod_strings,        only : integer_to_string
   use mod_communications, only : PAR_MAX
   use mod_par_affinity,   only : par_affinity

   implicit none

   integer(ip)                  :: i,j
   integer(ip),   pointer       :: app_id_arra(:)
   integer(ip)                  :: app_id,my_new_rank
   character(128)               :: app_name
   character(128), pointer      :: app_arra(:)
   character(128), pointer      :: app_arra_unique(:)
   integer(4)                   :: my_new_rank4,my_new_rank_wm4
   integer(ip)                  :: code_num, host_color, num_hosts
   integer(ip),   pointer       :: lcode(:)
   ! 
   ! Nullify pointers
   !
   nullify(app_arra)
   nullify(app_arra_unique)
   nullify(app_id_arra)
   !
   ! Get the names of the different Alya codes
   !
   call get_command_argument(COMMAND_ARGUMENT_COUNT(),app_name)

#ifndef MPI_OFF
   !
   ! Get your rank, and the world size. PAR_COMM is MPI_COMM_WORLD
   !
   call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_WORLD)
   call PAR_COMM_RANK_AND_SIZE(PAR_COMM_WORLD,PAR_MY_WORLD_RANK,PAR_WORLD_SIZE)

   if( PAR_WORLD_SIZE > 1 ) then
      !
      ! Rank 0 gathers the names, compares them and assigns app_id4 to each process
      !
      IPARALL = .true.
      if( PAR_MY_WORLD_RANK == 0 ) then
         allocate(app_id_arra(PAR_WORLD_SIZE))
         allocate(app_arra(PAR_WORLD_SIZE))
      end if

      call PAR_GATHER(app_name,app_arra,'IN THE WORLD')

      code_num = 1_ip ! default
      
      if( PAR_MY_WORLD_RANK == 0 ) then
         
         allocate(app_arra_unique(PAR_WORLD_SIZE))

         do i = 1_ip, PAR_WORLD_SIZE
            app_arra_unique(i) = app_arra(i)
         end do
         call string_sort(app_arra_unique, mcode,.false.)

         app_id_arra = -1_ip
         
         allocate(lcode(mcode))

         do i = 1_ip,mcode
#if defined COMMDOM && COMMDOM == 2
            code_num = 1_ip
#else
            code_num = read_code_number(app_arra_unique(i))
#endif

            if( code_num > mcode ) then
               print*,'PAR_CODE_SPLIT_WORLD: CODE '//integer_to_string(code_num)//' IN .DAT FILE EXCEED NUMBER OF ALYAs RUNNING'
               stop 1            
            end if

            lcode(code_num) = 1_ip

            do j = 1_ip, PAR_WORLD_SIZE
               if ( trim(app_arra_unique(i))==trim(app_arra(j)) ) then
                  app_id_arra(j) = code_num
               end if
            end do
         end do

         do j = 1_ip, PAR_WORLD_SIZE
            if( app_id_arra(j) < 0 ) then
               print*,'PAR_CODE_SPLIT_WORLD: ALYA '//trim(app_arra(j))//' TASK NUMBER '//trim(intost(j))//' WAS NOT ASSIGNED A CODE NUMBER'
               print*,'ALLOCATED CODE NUMBERS:'
               do i=1,mcode
                  print *,app_arra_unique(i),' - ',lcode(i)
               end do
               stop 1
            end if
         end do

         deallocate(app_arra_unique)
         deallocate(lcode)
      end if
      !
      ! Maximum number of codes
      !
      call PAR_BROADCAST(mcode,'IN THE WORLD')
      ! 
      ! Rank 0 scatters the app_id_arra4, and the split of COMM_WORLD is performed
      !
      call PAR_SCATTER(app_id_arra,app_id,'IN THE WORLD')
      !
      ! Split communicator
      !    
      call PAR_COMM_SPLIT(app_id,PAR_COMM_MY_CODE,my_new_rank,'IN THE WORLD')

      PAR_COMM_CURRENT  = PAR_COMM_MY_CODE 
      !
      ! My code communicator without master
      !
      my_new_rank4 = min(1_4,int(my_new_rank,4))
      call PAR_COMM_SPLIT(my_new_rank4,PAR_COMM_MY_CODE_WM,my_new_rank_wm4,'IN MY CODE')
      PAR_MY_CODE_RANK_WM = int(my_new_rank_wm4,ip)
      !
      ! Communicator in the host without the master
      !
      call PAR_DEFINE_COMMUNICATOR('IN THE UNIVERSE',PAR_COMM_UNIVERSE)
      call PAR_COMM_RANK_AND_SIZE(PAR_COMM_UNIVERSE,PAR_MY_UNIVERSE_RANK,PAR_UNIVERSE_SIZE)
      call par_affinity(host_num=host_color,num_hosts=num_hosts)
      if (PAR_MY_UNIVERSE_RANK >= 1) host_color = host_color + 1
      call PAR_COMM_SPLIT(host_color,PAR_COMM_SHARED_WM,PAR_SHARED_RANK_WM,PAR_COMM_UNIVERSE,PAR_MY_UNIVERSE_RANK)


      !
      ! Deallocate arrays 
      !
      if( PAR_MY_WORLD_RANK == 0 ) then
         deallocate(app_id_arra)
         deallocate(app_arra)
      end if

   else

      call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_MY_CODE)
     
      PAR_COMM_CURRENT = PAR_COMM_MY_CODE

   end if

#else
   
   PAR_WORLD_SIZE = 1

#endif
  
 contains

   !-----------------------------------------------------------------------
   !> 
   !> @author  Constantine Butakoff and houzeaux
   !> @date    2022-10-14
   !> @brief   Read code number
   !> @details Read the code number in data file
   !> 
   !-----------------------------------------------------------------------
   
   integer(ip) function read_code_number(name)
     
      use def_master,  only : lun_pdata
      use def_kintyp,  only : lg
      use mod_ecoute,  only : ecoute_set_read_unit
      use mod_ecoute,  only : ecoute
      use def_inpout,  only : getint,words
      use mod_iofile,  only : iofile_open_unit
      use mod_iofile,  only : iofile_close_unit
      implicit none 
      character(len=*), intent(in) :: name
      integer(ip)                  :: stat
      logical                      :: file_exists

      !Return 1 by default, if fails -- kill Alya
      !outside we expect that if CODE keyword was not found, we are running single code
      read_code_number = 1_ip
      
      INQUIRE(FILE=trim(name)//'.dat', EXIST=file_exists)
      if(.not. file_exists) then
         print *,'PAR_CODE_SPLIT_WORLD: FILE '//trim(name)//'.dat NOT FOUND'
         stop 1
      end if

      call iofile_open_unit(lun_pdata,trim(name)//'.dat','DATA FILE', IOSTAT=stat, crash_if_cannot_open=.false.)
      if(stat .ne. 0) then
         print *,'PAR_CODE_SPLIT_WORLD: FAILED TO OPEN '//trim(name)//'.dat. ERROR CODE '//trim(intost(stat))//' lun_pdata=',lun_pdata
         stop 1
      end if
      call ecoute_set_read_unit(lun_pdata)
      
      call ecoute('RRUDAT')
      if( words(1) /= 'RUNDA' ) call runend('RRUDAT: WRONG RUN_DATA CARD')
      do while( words(1) /= 'ENDRU' )
         if( words(1) == 'CODE ' ) read_code_number = getint('CODE ',1_ip,'#My code number')
         call ecoute('RRUDAT')
      end do
      
      call iofile_close_unit(lun_pdata)
      
    end function read_code_number
  
end subroutine par_code_split_world
