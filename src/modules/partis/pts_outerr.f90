!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_outerr()
  !------------------------------------------------------------------------
  !****f* master/pts_outerr
  ! NAME
  !    pts_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    Turnon
  !***
  !------------------------------------------------------------------------

  use def_master
  use def_kermod
  use def_domain
  use def_partis
  use mod_outfor,        only : outfor
  use mod_pts_injection, only : PTS_INJ_VELOC_GAUSSIAN
  use mod_pts_injection, only : PTS_INJ_VELOC_CONIC   
  use mod_pts_injection, only : PTS_INJ_VELOC_SPRAY   
  use mod_pts_injection, only : PTS_INJ_GEO_CIRCLE    
  use mod_pts_injection, only : PTS_INJ_GEO_ANNULUS   
  implicit none
  integer(ip)   :: ierro=0,iwarn=0,itype,iinj
!  character(20) :: messa

  ierro = 0
  iwarn = 0
  !
  ! Particle injection
  !
  do itype = 1,mtyla
     if( parttyp(itype) % kfl_exist /= 0 ) then
        if(  parttyp(itype) % kfl_grafo == 1 .or. &
             parttyp(itype) % kfl_buofo == 1 ) then
           if(  grnor == 0.0_rp .and. gravi(1) == 0.0_rp .and. &
                gravi(2) == 0.0_rp .and. gravi(3) == 0.0_rp ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,&
                   'GRAVITY IS MISSING FOR PARTICLE INJECTION')
           end if
        end if
     end if
     if( parttyp(itype) % kfl_brown /= 0 ) then
        if( kfl_prope == 0 .and. parttyp(itype) % diffu == 0.0_rp ) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,&
                'FOR BROWNIAN MOTION, VISCOSITY SHOULD BE DEFINED IN KERMOD')
        end if
     end if
  end do
  if( kfl_prope == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'PROPERTIES NOT DEFINED')
  end if
  !
  ! Thermodynamic model
  !
  if( kfl_thermo_pts == 0 ) then
     do itype = 1,mtyla
        if( parttyp(itype) % kfl_exist /= 0 ) then
           if( parttyp(itype) % kfl_therm /= 0 ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,&
                   'THERMODYNAMIC TRANSPORT SHOULD BE ON')
           end if
        end if
     end do
  end if
  !
  ! Element bin
  !
  if( kfl_usbin_pts /= 0 .and. kfl_element_bin == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'ACTIVATE BIN_ELEMENT OPTION IN KER.DAT FILE')
  end if
  !
  ! If there are sliup walls and normals are not allocated
  !
  if(IMASTER) then
     if( kfl_slip_wall_pts > 0 .or. kfl_bouncing_wall_pts > 0 ) then
        if( kfl_walln/=1 ) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,&
                'For slip walls in partis, add WALL_NORMALS to ker.dat to NUMERICAL_TREATMENT section')
        end if
     end if
  end if  
  !
  ! Check injectors
  !
  do iinj = 1,kfl_imax_injector_pts
     if( injection_pts(iinj) % kfl_geometry /= 0 ) then
        !
        ! Check type
        !
        if( injection_pts(iinj) % kfl_particle_type /= 0 ) then
           itype =  injection_pts(iinj) % kfl_particle_type
           if( parttyp(itype) % kfl_exist == 0 ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'NON-DEFINED PARTICLE TYPE IS INJECTED BY INJECTOR '//trim(intost(iinj)))
           end if
        end if
        !
        ! Check injection velocity strategies
        !
        if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_CONIC .or. &
            &injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_SPRAY .or. &
            &injection_pts(iinj) % kfl_veloc  == PTS_INJ_VELOC_GAUSSIAN) then
           !
           ! Conic velocity normal
           !
           if(  (injection_pts(iinj) % kfl_geometry /= PTS_INJ_GEO_CIRCLE) .and.&
              & (injection_pts(iinj) % kfl_geometry /= PTS_INJ_GEO_ANNULUS)  ) then
              call outfor(1_ip,momod(modul)%lun_outpu,'CONIC, SPRAY, AND GAUSSIAN VELOCITY MODELS ARE ONLY COMPATICLE WITH CIRLE OR ANNULUS INJECTOR IN INJECTOR '//trim(intost(iinj)))
           end if
        endif
     end if
  end do
  !
  ! Check if there is a wall distance criterion but distance to the wall does not exist
  !
  !if( dimin_pts > 0.0_rp .and. kfl_walld == 0 ) then
  !   ierro = ierro + 1
  !   call outfor(1_ip,momod(modul)%lun_outpu,&
  !        'DISTANCE TO THE WALL SHOULD BE COMPUTED BY KERMOD TO USE WALL-DISTANCE BASED DEPOSITION CRITERION')
  !end if
  !
  ! Stop
  !
  call errors(3_ip,ierro,iwarn,'NULL')

  !messa = intost(ierro)
  !if( ierro == 1 ) then
  !   call runend('PTS_OUTERR: '//trim(messa)//' ERROR HAS BEEN FOUND')
  !   call runend('PTS_OUTERR: '//trim(messa)//' ERROR HAS BEEN FOUND')
  !else if( ierro >= 2 ) then
  !   call runend('PTS_OUTERR: '//trim(messa)//' ERRORS HAVE BEEN FOUND')
  !end if

end subroutine pts_outerr
