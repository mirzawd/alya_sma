!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Communication_Toolbox
!> @{
!> @file    def_communications.f90
!> @author  houzeaux
!> @date    2020-05-13
!> @brief   Variables
!> @details Variables needed for communications
!-----------------------------------------------------------------------

module def_communications

  use def_kintyp_basic,      only : ip,rp,lg,r1p,i1p
  use def_master,            only : kfl_paral
  use def_master,            only : IPARALL,IMASTER,ISLAVE,INOTSLAVE,ISEQUEN
  use def_master,            only : INOTMASTER
  use def_kintyp_comm,       only : comm_data_par
  use def_kintyp_comm,       only : comm_data_par_basic
  use def_master,            only : current_code,current_zone,current_subd
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_memory_basic,      only : memory_resize
  use mod_parall,            only : par_code_zone_subd_to_color
  use mod_parall,            only : PAR_COMM_COLOR
  use mod_parall,            only : PAR_COMM_MY_CODE
  use mod_parall,            only : PAR_COMM_MY_CODE_WM
  use mod_parall,            only : PAR_COMM_SHARED_WM
  use mod_parall,            only : PAR_COMM_COLOR_ARRAY
  use mod_parall,            only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,            only : PAR_MY_CODE_RANK
  use mod_parall,            only : PAR_MY_CODE_RANK_WM
  use mod_parall,            only : PAR_SHARED_RANK_WM
  use mod_parall,            only : PAR_MY_WORLD_RANK
  use mod_parall,            only : PAR_COMM_COLOR_PERM
  use mod_parall,            only : PAR_COMM_NULL
  use mod_parall,            only : PAR_COMM_WORLD
  use mod_parall,            only : PAR_COMM_CURRENT
  use mod_parall,            only : PAR_INTEGER
  use mod_parall,            only : PAR_REAL
  use mod_parall,            only : PAR_COMM_SFC
  use mod_parall,            only : PAR_COMM_SFC_WM
  use mod_parall,            only : PAR_COMM_MPIO
  use mod_parall,            only : PAR_COMM_MPIO_WM
  use mod_parall,            only : commd
  use mod_parall,            only : color_target
  use mod_parall,            only : color_source
  use mod_parall,            only : par_memor
  use mod_parall,            only : PAR_CODE_SIZE
  use mod_optional_argument, only : optional_argument
  use mod_std  
  use def_mpi

end module def_communications
