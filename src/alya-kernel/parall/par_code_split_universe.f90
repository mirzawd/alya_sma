!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!>
!> @addtogroup Parall
!> @{
!> @author  J.C. Cajas
!> @date    13/05/2014
!> @brief   Split MPI_COMM_WORLD in application communicators
!> @details Split MPI_COMM_WORLD in application communicators, based on the
!>          name of the application launched
!>
!>                        MPI_COMM_WORLD = PAR_COMM_UNIVERSE                                   |
!>                                 (THE UNIVERSE)                                              |
!>                                       ||                                                    |
!>                                       \/                                                    |
!>                                                                                             |
!>                                 MPI_COMM_SPLIT                                              |
!>                                                                                             |
!>                                //            \\                                             |           PERFORMED IN
!>                               //              \\                                            |    par_code_split_universe.f90
!>                                                                                             |
!>                   PAR_COMM_WORLD             REST OF THE WORLD                              |
!>              (ALYA's COMMUNICATOR)        (OTHER CODES COMMUNICATOR)                        |
!>                        ||                                                                   |
!>                        \/                                                                   |
!> ---------------------------------------------------------------------------------------------
!>                  MPI_COMM_SPLIT                                                             |
!>                                                                                             |
!>                        ||                                               ||                  |
!>                        \/                                               \/                  |
!>                                                                                             |           PERFORMED IN
!> PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD       PAR_COMM_COLOR_ARRAY(:) % PAR_COMM_WORLD   |      par_code_split_world.f90
!>                  = PAR_COMM_MY_CODE                                                         |   par_zone_communication_arrays.f90
!>                 (par_sencom.f90)                          (par_color_communicators.f90)     |
!>                                                                                             |
!>                                                                                             |
!> @}                                                                                          |
!----------------------------------------------------------------------------------------------+

subroutine par_code_split_universe()
  
   use def_kintyp,         only : ip
   use def_master,         only : IPARALL
   use def_master,         only : application_names
   use mod_parall,         only : PAR_COMM_WORLD
   use mod_parall,         only : PAR_COMM_UNIVERSE
   use mod_parall,         only : PAR_UNIVERSE_SIZE
   use mod_parall,         only : PAR_MY_UNIVERSE_RANK
   use mod_parall,         only : PAR_COMM_CURRENT
   use mod_parall,         only : mapps
   use mod_communications, only : PAR_INIT
   use mod_communications, only : PAR_COMM_RANK_AND_SIZE
   use mod_communications, only : PAR_GATHER
   use mod_communications, only : PAR_COMM_SPLIT
   use mod_communications, only : PAR_SCATTER
   use mod_communications, only : PAR_BROADCAST
   use mod_communications, only : PAR_DEFINE_COMMUNICATOR
   use def_mpi
#if defined COMMDOM && COMMDOM == 2
    use mod_plepp_pdn_contact, only : plepp_pdn_init
#endif
   use mod_strings,        only : string_sort
  
   implicit none

   integer(ip)             :: i,ii,j
   integer(ip)             :: path_f, ext_i  ! Variables to locate the end of the path and the beginning of the extension in app_type
   integer(ip)             :: app_type_id,my_alya_rank
   character(128)          :: app_type
   character(128), pointer :: app_type_arra(:)
   character(128), pointer :: app_name(:)
   integer(ip),    pointer :: app_type_id_arra(:)

#ifdef _WIN32 
      character(1), parameter :: path_sep = "\"
#else
      character(1), parameter :: path_sep = "/"
#endif

   !
   ! Initialize variables and pointers
   !
   call par_initialize_nullify()
   !
   ! Initialize MPI
   !
   call PAR_INIT()
   !
   ! Nullify pointers
   !
   nullify(app_type_id_arra)
   nullify(app_type_arra)
   nullify(app_name)
   !
   ! Get the types of the different applications launched
   ! Remove the path of the executable and the extension
   !
   call GET_COMMAND_ARGUMENT(0_4,app_type)

   path_f = 0_ip
   ext_i  = len_trim(app_type)+1_ip
   do ii = 1_ip, len_trim(app_type)      
      if( app_type(ii:ii) == path_sep )then
         path_f = ii
      end if
   end do
   do ii = len_trim(app_type),1_ip,-1_ip
      if( app_type(ii:ii) == '.' )then
         ext_i = ii
         exit
      end if
   end do
   app_type = trim(app_type(path_f+1_ip:ext_i-1_ip))

#ifndef MPI_OFF
   !
   ! Get your rank, and the world size. PAR_COMM is MPI_COMM_WORLD
   !
   call PAR_DEFINE_COMMUNICATOR('IN THE UNIVERSE',PAR_COMM_UNIVERSE)
   call PAR_COMM_RANK_AND_SIZE(PAR_COMM_UNIVERSE,PAR_MY_UNIVERSE_RANK,PAR_UNIVERSE_SIZE)

#if defined COMMDOM && COMMDOM == 2
   call plepp_pdn_init()
#else
   if( PAR_UNIVERSE_SIZE > 1 ) then
      !
      ! Rank 0 gathers the types, compares them and assigns type_id4 to each process
      !
      IPARALL = .true.
      if( PAR_MY_UNIVERSE_RANK == 0 ) then
         allocate(app_type_arra(PAR_UNIVERSE_SIZE))
         allocate(app_type_id_arra(PAR_UNIVERSE_SIZE))
         allocate(app_name(PAR_UNIVERSE_SIZE))
      end if

      call PAR_GATHER(app_type,app_type_arra,'IN THE UNIVERSE')
      
      if( PAR_MY_UNIVERSE_RANK == 0 ) then ! If I am de MASTER OF THE UNIVERSE!!!

         do i = 1_ip, PAR_UNIVERSE_SIZE
            app_name(i) = app_type_arra(i)
         end do
         call string_sort(app_name, mapps,.false.)

         app_type_id_arra(:) = -1
         do i = 1_ip, mapps
            do j = 1_ip, PAR_UNIVERSE_SIZE
               if ( trim(app_name(i))==trim(app_type_arra(j)) ) then
                  app_type_id_arra(j) = i
               end if
            end do
         end do
      end if

      if( PAR_MY_UNIVERSE_RANK == 0 ) then
         allocate(application_names(mapps))
         do i = 1,mapps
            application_names(i) = trim(app_name(i))
         end do
      end if
      !
      ! Maximum number of codes
      !
      call PAR_BROADCAST(mapps,'IN THE UNIVERSE')
      !
      ! Rank 0 scatters the type_id_arra4, and the split of COMM_WORLD in applications is performed
      !
      call PAR_SCATTER(app_type_id_arra,app_type_id,'IN THE UNIVERSE')

      if( mapps > 1_ip ) then
         !
         ! Split communicator. This defines the communicator PAR_COMM_WORLD (the Alya world)
         !
         call PAR_COMM_SPLIT(app_type_id,PAR_COMM_WORLD,my_alya_rank,'IN THE UNIVERSE')        
         PAR_COMM_CURRENT  = PAR_COMM_WORLD

      else

         PAR_COMM_WORLD    = MPI_COMM_WORLD
         PAR_COMM_CURRENT  = PAR_COMM_WORLD

      end if
      !
      ! Deallocate arrays
      !
      if( PAR_MY_UNIVERSE_RANK == 0 ) then

         deallocate(app_type_id_arra)
         deallocate(app_type_arra)
         deallocate(app_name)

      end if

   else
      PAR_COMM_WORLD    = MPI_COMM_WORLD
      PAR_COMM_CURRENT  = PAR_COMM_WORLD
   end if
#endif

#endif

end subroutine par_code_split_universe
