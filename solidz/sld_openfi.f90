!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_openfi.f90
!> @author  Gerard Guillamet
!> @date    July, 2018
!> @brief   Open file names
!> @details
!>
!>          This subroutine gets ALL the file names and open them to be
!>          used by the module in two possible ways:\n
!>
!>          (a) Recalling them from the environment, when Alya is
!>              launched encapsulated in a shell script, or\n
!>          (b) Composing the names out of the problem name which is
!>              given as argument when the binary file Alya is
!>              launched "naked".\n
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_openfi(itask)

  use def_kintyp,  only : ip, rp
  use def_master,  only : INOTSLAVE, modul, namda, exmod
  use def_master,  only : kfl_naked, kfl_rstar
  use mod_iofile,  only : iofile
  use def_solidz,  only : kfl_psmat_sld, lun_psmat_sld
  use def_solidz,  only : lun_react_sld
  use def_solidz,  only : lun_carcy_res_sld, lun_carcy_cvg_sld
  use def_solidz,  only : lun_sysnet_heart_res_sld, lun_sysnet_system_res_sld
  use mod_sld_cardiac_cycle,     only : kfl_cardiac_cycle
  use mod_sld_post_reaction,     only : kfl_preac_sld
#if defined COMMDOM && COMMDOM == 2
  use def_solidz,                only : kfl_rigid_sld, kfl_conta_sld, kfl_contf_sld
  use mod_sld_pdn_contact_plepp, only : kfl_pdnco_sld, lun_pdnco_sld, kfl_pdncf_sld
#endif
#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet,            only : kfl_sysnet
#endif

  implicit none

  integer(ip),   intent(in) :: itask
  character(150)            :: fil_psmat_sld
  character(150)            :: fil_react_sld
  character(150)            :: fil_pdnco_sld
  character(150)            :: fil_carcy_res_sld, fil_sysnet_heart_res_sld
  character(150)            :: fil_carcy_cvg_sld, fil_sysnet_system_res_sld
  character(7)              :: statu
  character(11)             :: forma
  character(6)              :: posit

  if ( INOTSLAVE ) then
     !
     ! Define unit opening option if this is a restart run
     !
     if(kfl_rstar == 2_ip) then
        statu='old'
        forma='formatted'
        posit='append'
     else
        statu='unknown'
        forma='formatted'
        posit='asis'
     end if

     select case (itask)

     case (1_ip)
        !
        ! Open files needed occasionally
        !
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR114',fil_psmat_sld)
           call GET_ENVIRONMENT_VARIABLE('FOR115',fil_carcy_res_sld)
           call GET_ENVIRONMENT_VARIABLE('FOR116',fil_carcy_cvg_sld)
           call GET_ENVIRONMENT_VARIABLE('FOR117',fil_react_sld)
           call GET_ENVIRONMENT_VARIABLE('FOR118',fil_pdnco_sld)
        else
           fil_psmat_sld = adjustl(trim(namda))//'-matrix.'//exmod(modul)//'.ps'
           fil_carcy_res_sld = adjustl(trim(namda))//'-cardiac-cycle.'//exmod(modul)//'.res'
           fil_carcy_cvg_sld = adjustl(trim(namda))//'-cardiac-cycle.'//exmod(modul)//'.cvg'
           fil_react_sld = adjustl(trim(namda))//'-reaction.'//exmod(modul)//'.res'
           fil_sysnet_heart_res_sld = adjustl(trim(namda))//'-sysnetHeart.'//exmod(modul)//'.res'
           fil_sysnet_system_res_sld = adjustl(trim(namda))//'-sysnetSystem.'//exmod(modul)//'.res'
           fil_pdnco_sld = adjustl(trim(namda))//'-contact.'//exmod(modul)//'.cvg'
        end if
        !
        ! Matrix profile
        !
        if( kfl_psmat_sld > 0 ) then
           call iofile(0_ip,lun_psmat_sld,fil_psmat_sld,'SOLIDZ MATRIX',statu,forma,posit)
        end if
        !
        ! Cardiac cycle
        !
        if( kfl_cardiac_cycle ) then
           call iofile(0_ip,lun_carcy_res_sld,fil_carcy_res_sld,'SOLIDZ CARDIAC CYCLE',statu,forma,posit)
           call iofile(0_ip,lun_carcy_cvg_sld,fil_carcy_cvg_sld,'SOLIDZ CARDIAC CYCLE',statu,forma,posit)
        end if
        !
        ! Displacements and reactions on set of nodes
        !
        if( kfl_preac_sld ) then
           call iofile(0_ip,lun_react_sld,fil_react_sld,'SOLIDZ REACTIONS ON NODES DEFINED BY THE USER',statu,forma,posit)
        end if
        !
        ! PDN-contact convergence
        !
#if defined COMMDOM && COMMDOM == 2
        if( (kfl_pdnco_sld > 0 .or. kfl_conta_sld == 3_ip) .and. &
            (kfl_pdncf_sld .or. kfl_contf_sld == 1) .and. kfl_rigid_sld == 0_ip ) then
           call iofile(0_ip,lun_pdnco_sld,fil_pdnco_sld,'SOLIDZ PDN CONTACT CONVERGENCE',statu,forma,posit)
        end if
#endif
#ifndef PROPER_ELEM_PRIVATE_OFF
            !
            ! Sysnet
            !
            if( kfl_sysnet ) then
                call iofile(0_ip,lun_sysnet_heart_res_sld,fil_sysnet_heart_res_sld,'SOLIDZ CARDIAC CYCLE',statu,forma,posit)
                call iofile(0_ip,lun_sysnet_system_res_sld,fil_sysnet_system_res_sld,'SOLIDZ CARDIAC CYCLE',statu,forma,posit)
            end if
#endif

     case (2_ip)
        !
        ! Close files
        !
        call iofile(2_ip,lun_psmat_sld,' ','SOLIDZ MATRIX PROFILE')
        call iofile(2_ip,lun_react_sld,' ','SOLIDZ REACTIONS ON NODES DEFINED BY THE USER')
#if defined COMMDOM && COMMDOM == 2
        call iofile(2_ip,lun_pdnco_sld,' ','SOLIDZ CONVERGENCE PDN-CONTACT')
#endif

     end select

  end if

end subroutine sld_openfi

