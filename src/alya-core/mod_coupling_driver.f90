!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling functions
!> @file    mod_coupling_driver.f90
!> @author  Guillaume Houzeaux
!> @date    11/06/2014
!> @brief   Driver for coupling
!> @details Driver for coupling
!> @{
!------------------------------------------------------------------------

module mod_coupling_driver
  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : intost
  use def_master,         only : modul,momod
  use def_master,         only : lmord
  use def_master,         only : mmodu
  use def_master,         only : ittim
  use def_master,         only : iblok
  use def_master,         only : ittim
  use def_master,         only : IMASTER,ITASK_INIUNK
  use def_master,         only : ID_NASTIN
  use def_master,         only : ID_TEMPER
  use def_master,         only : ID_SOLIDZ
  use def_master,         only : ID_ALEFOR
  use def_master,         only : ID_PARTIS
  use def_master,         only : ID_CHEMIC
  use def_master,         only : ID_KERMOD
  use def_master,         only : ID_TURBUL
  use def_master,         only : ID_EXMEDI
  use def_master,         only : ID_GUSANO
  use def_master,         only : ID_KERNEL
  use def_master,         only : ITASK_AFTER
  use def_master,         only : namod
  use mod_parall,         only : I_AM_IN_COLOR
  use def_coupli,         only : coupling_type
  use def_coupli,         only : BETWEEN_ZONES
  use def_coupli,         only : mcoup
  use def_coupli,         only : coupling_driver_couplings
  use def_coupli,         only : coupling_driver_number_couplings
  use mod_ker_timeline,   only : ker_timeline  
  use mod_messages,       only : messages_live
  use mod_alya2dlb,       only : alya2dlb_DLB_Barrier
  use mod_moduls_conf,    only : moduls_in_block
  use def_master,         only : times
  use mod_module_interface
#if defined COMMDOM && COMMDOM == 2
  use mod_plepp_pdn_contact, only : CNT_CPLNG
  use mod_plepp_pdn_contact, only : plepp_pdn_driver_sendrecv
#endif
  
  implicit none

contains

  subroutine COU_DRIVER(current_when,current_task)
    integer(ip),  intent(in)  :: current_when
    integer(ip),  intent(in)  :: current_task
#ifndef COMMDOM
    integer(ip)               :: icoup,kcoup
    integer(ip)               :: ITASK_COUPL
    integer(ip)               :: module_source
    integer(ip)               :: module_target
    integer(ip)               :: color_source
    integer(ip)               :: color_target
    logical(lg)               :: i_compute_and_send
    logical(lg)               :: i_recv_and_assemble
#endif
#ifdef ALYA_DLB_BARRIER
    integer(ip)               :: ierror
#endif
    !
    ! Loop over couplings
    !
#ifndef COMMDOM

    call momod(0) % times(11) % ini()
    
    do kcoup = 1,coupling_driver_number_couplings(iblok)
       icoup = coupling_driver_couplings(kcoup,iblok)

       if( coupling_type(icoup) % kind == BETWEEN_ZONES ) then

          module_source = coupling_type(icoup) % module_source
          module_target = coupling_type(icoup) % module_target
          color_source  = coupling_type(icoup) % color_source
          color_target  = coupling_type(icoup) % color_target
          ITASK_COUPL   = icoup + 1000
          !
          ! Should I stay or should I go
          !
          ! Only activate coupling if the current time step is a multiple of the coupling frequency
          !
          i_compute_and_send  = .false.
          if(    current_task == coupling_type(icoup) % task_compute_and_send .and. &
               & current_when == coupling_type(icoup) % when_compute_and_send .and. &
               & I_AM_IN_COLOR(color_source).and. modul == module_source      .and. &
               & mod( ittim,coupling_type(icoup) % frequ_send ) == 0_ip )   then
             i_compute_and_send  = moduls_in_block(iblok,modul) 
          end if
          
          i_recv_and_assemble = .false.
          if(    current_task == coupling_type(icoup) % task_recv_and_assemble .and. &
               & current_when == coupling_type(icoup) % when_recv_and_assemble .and. &
               & I_AM_IN_COLOR(color_target).and. modul == module_target       .and. &
               & mod( ittim,coupling_type(icoup) % frequ_recv ) == 0_ip)     then
             i_recv_and_assemble = moduls_in_block(iblok,modul) 
          end if
          !
          ! Call corresponding module (should call the plugin of the module_source/target and that's it)
          !
          if( module_source == modul .or. module_target == modul ) then
             if(  ( i_compute_and_send  .and. .not. i_recv_and_assemble ) .or. &
                & ( i_recv_and_assemble .and. .not. i_compute_and_send  ) ) then
                !
                ! I am source or target
                !
                if( i_compute_and_send ) then
                   call messages_live('<- SEND '//trim(coupling_type(icoup) % variable)//' AS A SOURCE FOR COUPLING '//trim(intost(icoup)))
                   call ker_timeline(module_source,'INI_COUPLING',icoup,module_source)
                else if( i_recv_and_assemble ) then
                   call messages_live('-> RECV '//trim(coupling_type(icoup) % variable)//' AS A TARGET FOR COUPLING '//trim(intost(icoup)))
                end if
                !
                ! DLB Barrier, to calm down MPI...
                !
#ifdef ALYA_DLB_BARRIER
                if( modul > 0 ) then
                   if( IMASTER ) call messages_live(namod(modul)//': BEFORE DLB BARRIER')
                   ierror = alya2dlb_DLB_Barrier()
                   if( IMASTER ) call messages_live(namod(modul)//': AFTER DLB BARRIER')
                end if
#endif
#ifdef ALYA_EXTRAE
                if( modul > 0 ) call extrae_eventandcounters(900,99_8) 
#endif

                select case ( modul )

                case ( ID_NASTIN )
                   call Nastin(ITASK_COUPL)
                case ( ID_TEMPER )
                   call Temper(ITASK_COUPL)
                case ( ID_SOLIDZ )
                   call Solidz(ITASK_COUPL)
                case ( ID_TURBUL )
                   call Turbul(ITASK_COUPL)
                case ( ID_ALEFOR )
                   call Alefor(ITASK_COUPL)
                case ( ID_PARTIS )
                   call Partis(ITASK_COUPL)
                case ( ID_CHEMIC )
                   call Chemic(ITASK_COUPL)
                case ( ID_KERMOD )
                   call Kermod(ITASK_COUPL)
                case ( ID_EXMEDI )
                   call Exmedi(ITASK_COUPL)
                case ( ID_GUSANO )
                   call Gusano(ITASK_COUPL)
                case default
                   call runend('COU_DRIVER: MODULE NOT CODED')

                end select

#ifdef ALYA_EXTRAE
                if( modul > 0 ) call extrae_eventandcounters(900,0_8) 
#endif
                
             else if( i_compute_and_send .and. i_recv_and_assemble ) then

                call runend('COU_DRIVER: WE ARE IN TROUBLE')

             end if

          end if

          if( i_compute_and_send ) then
             call ker_timeline(module_source,'END_COUPLING',icoup,module_target)
          end if

       end if

    end do
    
#else
    
#if COMMDOM == 2
    call plepp_pdn_driver_sendrecv( CNT_CPLNG, current_when, current_task )
#endif

#endif
    call momod(0) % times(11) % add()

  end subroutine COU_DRIVER

end module mod_coupling_driver
!> @}
!-----------------------------------------------------------------------
