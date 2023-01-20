!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine outerr(itask)
  !------------------------------------------------------------------------
  !****f* master/outerr
  ! NAME 
  !    outerr
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
  use mod_parall,   only : PAR_PARALLEL_PARTITION
  use mod_parall,   only : PAR_METIS4
  use def_parall,   only : kfl_partition_par
  use def_parall,   only : kfl_parseq_par
  use mod_outfor,   only : outfor
  use def_AMR,      only : kfl_amr
  use mod_messages, only : messages_live
  use mod_live_info_config, only : live_info_config
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ierro=0,iwarn=0,imodu,iorde
  integer(ip)             :: korde(mmodu),kblok
  character(20)           :: messa

  ierro = 0
  iwarn = 0

  select case ( itask )

  case ( 0_ip )
     !
     ! Velocity function
     !
     if( kfl_vefun /= 0 .and. kfl_modul(ID_NASTIN) + kfl_modul(ID_NASTAL) > 0 ) then
        ierro = ierro + 1
        call outfor(1_ip,live_info_config%lun_livei,&
             'CANNOT SOLVE NASTIN OR NASTAL AND USE A VELOCITY FUNCTION')        
     end if
     !
     ! METIS
     !
#ifdef METISI8
#ifdef V5METIS
     ierro = ierro + 1
     call outfor(1_ip,live_info_config%lun_livei,&
          'METISI8 AND V5METIS NOT COMPATIBLE')        
#endif
#endif
     !
     ! METIS in parallel
     !
     if( kfl_parseq_par == PAR_PARALLEL_PARTITION .and. kfl_partition_par == PAR_METIS4 ) then
        ierro = ierro + 1
        call outfor(1_ip,live_info_config%lun_livei,&
             'METISI4 SHOULD BE USED IN SEQUENTIAL EXECUTION MODE, CHECK THE .dat file')        
     end if


  case ( 1_ip )

     if( INOTSLAVE ) then
        !
        ! Check if a problem is put to on that it is solved
        !
        korde=0
        do kblok=1,nblok
           do imodu=1,mmodu
              iorde=lmord(imodu,kblok)
              if(iorde>0) then
                 korde(iorde)=korde(iorde)+1
              end if
           end do
        end do
        do imodu=1,mmodu-1
           if(kfl_modul(imodu)/=0.and.korde(imodu)==0) then
              ierro=ierro+1
              messa=intost(imodu)
              call outfor(1_ip,live_info_config%lun_livei,&
                   'MODULE '//trim(messa)//' IS SOLVED BUT DOES NOT APPEAR'//&
                   ' IN BLOCK DEFINITION')
           end if
        end do
        !
        ! Witness points
        !
        if(nwitn/=0.and.kfl_elses==0) then
           ierro=ierro+1
           nwitn=0
           call outfor(1_ip,live_info_config%lun_livei,&
                'DEFINE ELSEST FIELD TO OUTPUT WITNESS POINTS')
        end if
        
        !
        ! Periodicity does not work with SFC 
        !
        if( nperi > 0_ip .and. kfl_partition_par /= PAR_METIS4 ) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,&
                'FOR PERIODIC CASES YOU NEED TO USE METHOD: METIS IN THE .DAT')
        end if
     end if
     !
     ! AMR and mesh multiplication
     !
     if( kfl_amr /= 0 .and. ndivi > 0 ) then
        ierro = ierro + 1
        call messages_live('ERROR: AMR AND MESH MULTIPLICATION NOT COMPATIBLE')
     end if
     
  end select
  !
  ! Stop
  !
  messa = intost(ierro)
  if( ierro == 1 ) then
     call runend('OUTERR: '//trim(messa)//' ERROR HAS BEEN FOUND')
  else if( ierro >= 2 ) then
     call runend('OUTERR: '//trim(messa)//' ERRORS HAVE BEEN FOUND')
  end if

end subroutine outerr
