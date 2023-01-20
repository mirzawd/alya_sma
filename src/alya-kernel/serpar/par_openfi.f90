!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_openfi(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_openfi
  ! NAME
  !    par_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names to be used by Parall 
  !    service of Alya in two
  !    possible ways and them open them:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked". 
  !
  !    The files to be opened are:
  !
  !    MASTER
  !    lun_domai_par ... Domain partition
  !    lun_outpu_par ... Output general information
  !    SLAVES
  !    lun_memor ....... Memory information
  !
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_parall
  use mod_iofile
  use mod_parall
  use mod_memory, only : lun_memor
  use mod_par_virfil
  use mod_outfor, only : outfor
  use mod_run_config, only : run_config
  use mod_memory_config, only : memory_config
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iunit
  character(150)          :: fil_domai_par,fil_outpu_par
  character(150)          :: fil_trace_par,fil_memor_par
  character(150)          :: fil_postp_par,fil_conve_par,cfile
  character(150)          :: fil_matri_msh,fil_matri_res 
  character(150)          :: fil_repar_par
  character(20)           :: cdum1,cdum2

  select case (itask)

  case(0_ip)

     !-------------------------------------------------------------------
     !
     ! Preprocess and Restart
     !
     !-------------------------------------------------------------------
     
     if( ( PART_AND_WRITE() .or. READ_AND_RUN() ) .and. IPARALL ) then
        if( kfl_naked == 0 ) then
           call GET_ENVIRONMENT_VARIABLE('FOR5516',fil_rstar_par)
        else if( kfl_naked == 1 ) then
           fil_rstar_par = adjustl(trim(namda))//'.par.rst'
        end if

        if( PART_AND_WRITE() .and. nproc_par > 1 ) then
           !
           ! Parallel preprocess
           !
           call par_filnam(1_ip,kfl_paral,fil_rstar_par,cfile)
           iunit = lun_aonlp_par 
           cdum1 = intost(kfl_paral)
           cdum2 = 'PARALL RESTART '//trim(cdum1)
           call iofile(zero,iunit,trim(cfile),trim(cdum2),'replace')
           close(iunit)
           if( kfl_filio_par == 0 .and. kfl_ascii_par == 0 ) then
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','unformatted','append')
           else if( kfl_filio_par == 0 ) then
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','formatted','append')                 
           end if
           if( kfl_virfi_par == 1 ) call par_inibuf(-1_ip)
           
        else if( PART_AND_WRITE() .and. IMASTER ) then
           !
           ! Sequential Preprocess
           !
           do kfl_desti_par = 0,npart_par
              call par_filnam(1_ip,kfl_desti_par,fil_rstar_par,cfile)
              iunit = lun_aonlp_par + kfl_desti_par
              cdum1 = intost(kfl_desti_par)
              cdum2 = 'PARALL RESTART '//trim(cdum1)
              call iofile(zero,iunit,trim(cfile),trim(cdum2),'replace')
              close(iunit)
              if( kfl_filio_par == 0 .and. kfl_ascii_par == 0 ) then
                 call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','unformatted','append')
              else if( kfl_filio_par == 0 ) then
                 call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','formatted','append')                 
              end if
           end do
           if( kfl_virfi_par == 1 ) call par_inibuf(-1_ip)

        else if( READ_AND_RUN() .and. IPARALL ) then
           !
           ! Restart
           ! 
           call par_filnam(1_ip,kfl_paral,fil_rstar_par,cfile)
           iunit = lun_aonlp_par !+ kfl_paral
           
           if( kfl_ascii_par == 0 ) then
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','unformatted')
           else
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old')
           end if 

        end if
     end if

  case(1_ip)

     !-------------------------------------------------------------------
     !
     ! Get file names:
     !
     ! If kfl_naked=0 -->> encapsulated, then get names from the environment (DEFAULT value)
     ! If kfl_naked=1 -->> naked, then compose the names
     !    
     !-------------------------------------------------------------------

     if(kfl_naked==0) then
        call GET_ENVIRONMENT_VARIABLE('FOR5502',fil_outpu_par) 
        call GET_ENVIRONMENT_VARIABLE('FOR5504',fil_trace_par) 
        call GET_ENVIRONMENT_VARIABLE('FOR5505',fil_memor_par)
        call GET_ENVIRONMENT_VARIABLE('FOR5506',fil_conve_par)
        call GET_ENVIRONMENT_VARIABLE('FOR5507',fil_parti_msh)
        call GET_ENVIRONMENT_VARIABLE('FOR5508',fil_parti_res)
        call GET_ENVIRONMENT_VARIABLE('FOR5509',fil_matri_msh)
        call GET_ENVIRONMENT_VARIABLE('FOR5510',fil_matri_res)
        call GET_ENVIRONMENT_VARIABLE('FOR5511',fil_conne_par)
        call GET_ENVIRONMENT_VARIABLE('FOR5512',fil_repar_par)
        call GET_ENVIRONMENT_VARIABLE('FOR5515',fil_postp_par)
        call GET_ENVIRONMENT_VARIABLE('FOR5516',fil_rstar_par)
     else if(kfl_naked==1) then
        fil_outpu_par  = adjustl(trim(namda))//'.par.log'
        fil_trace_par  = adjustl(trim(namda))//'.par.tra'
        fil_memor_par  = adjustl(trim(namda))//'.par.mem'
        fil_conve_par  = adjustl(trim(namda))//'.par.cvg'
        fil_parti_msh  = adjustl(trim(namda))//'-partition.par.post.msh'
        fil_parti_res  = adjustl(trim(namda))//'-partition.par.post.res'
        fil_matri_msh  = adjustl(trim(namda))//'-matrix.par.post.msh'
        fil_matri_res  = adjustl(trim(namda))//'-matrix.par.post.res'
        fil_conne_par  = adjustl(trim(namda))//'-connectivity.par.res'
        fil_repar_par  = adjustl(trim(namda))//'-repartitioning.par.res'
        fil_postp_par  = adjustl(trim(namda))//'.par.post.res'
        fil_rstar_par  = adjustl(trim(namda))//'.par.rst'
     end if
     !
     ! Open master files
     !
     if( IMASTER ) then
        call iofile(zero,lun_outpu_par,fil_outpu_par,'PARALL OUTPUT')
        call iofile(zero,lun_conve_par,fil_conve_par,'PARALL COMMUNICATION TIME')
        if( kfl_matri_par == 1 ) then
           call iofile(zero,lun_matri_msh,fil_matri_msh,'PARALL MATRIX MESH')
           call iofile(zero,lun_matri_res,fil_matri_res,'PARALL MATRIX RESULT')           
        end if
     end if
     !
     ! Output (log) file
     !
     if( ISLAVE .and. .not. PART_AND_WRITE() ) then
        memory_config%output = .false.
     end if
     !
     ! Postprocess is performed by slaves
     !
     if( ISLAVE .and. .not. PART_AND_WRITE() ) then
        if( kfl_postp_par == 0 ) then
           fil_postp_par = trim(fil_postp_par)//trim(intost(kfl_paral))
           call iofile(zero,lun_postp,fil_postp_par,'POST-PROCESS')
           write(lun_postp,'(a)')'GiD Post Results File 1.0'
           write(lun_postp,'(a)')' '
        end if
     end if
     !
     ! repartitioning
     !
     if( kfl_repart_par /= 0 .and. IMASTER ) then
        if( kfl_rstar == 2 ) then
           call iofile(zero,lun_repar_par,trim(fil_repar_par),'REPARTITIONING','old','formatted','append')
        else
           call iofile(zero,lun_repar_par,trim(fil_repar_par),'REPARTITIONING')
        end if
     end if
     
  case(2_ip)
     !
     ! Open domain file
     !     
     if( kfl_naked == 0 ) then
        call GET_ENVIRONMENT_VARIABLE('FOR5503',fil_domai_par)
     else if( kfl_naked == 1 ) then
        fil_domai_par = adjustl(trim(namda))//'.par.post.msh'
     end if
     if( ISLAVE .and. kfl_postp_par == 0 ) then
        if( kfl_oumes(1) == 1 ) then
           fil_domai_par = trim(fil_domai_par)//trim(intost(kfl_paral))
           call iofile(zero,lun_domai_par,fil_domai_par,'PARALL DOMAIN') 
        end if
     end if

  case(4_ip)

     !-------------------------------------------------------------------
     !
     ! Close files
     !
     !-------------------------------------------------------------------

     if( IMASTER ) then
        call iofile(two,lun_outpu_par,' ','PARALL OUTPUT')        
     end if

  case(5_ip)

     !-------------------------------------------------------------------
     !
     ! Compose file name according to master/slave number
     !
     !-------------------------------------------------------------------

     if( IPARALL ) then
        parch=trim(parch)//'.par'//trim(intost(kfl_paral))
     end if
     
  end select

end subroutine par_openfi
