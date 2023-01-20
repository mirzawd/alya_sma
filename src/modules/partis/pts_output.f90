!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @name    Starts an iteration
!> @file    pts_output.f90
!> @author  Guillaume Houzeaux
!> @date    28/01/2013
!> @brief   Output of particles
!> @details Postprocess on mesh (particle density) and output
!> @}
!------------------------------------------------------------------------

subroutine pts_output()
  use def_parame
  use def_master
  use def_kermod
  use def_partis
  use mod_memory
  use mod_ker_timeline,        only : ker_timeline
  use mod_communications,      only : PAR_GATHER,PAR_GATHERV
  use mod_pts_parallelization, only : pts_parallelization_pack
  use mod_pts_parallelization, only : pts_parallelization_unpack
  use mod_messages,            only : livinf
  use mod_outfor,              only : outfor
  use mod_iofile,              only : iofile_flush_unit
  use mod_output_postprocess,  only : output_postprocess_variables
  use mod_result_io,           only : pts_result_io
  
  implicit none

  external                 :: pts_outvar
  integer(ip)              :: ivari,ivarp,ilagr,ipars,ipart,ivard
  integer(ip)              :: kfl_ifpos,kfl_ifdep,imaxi,imini,nvart_pts
  integer(4)               :: nvarp_pts4,nvard_pts4,nvart_pts4
  integer(ip), save        :: ittim_last=-1
  real(rp)                 :: cutim_pts
  real(rp),    allocatable :: depos_surface(:)
  integer(4)               :: num_particles
  integer(4),  pointer     :: num_particles_gat(:)
  real(rp),    pointer     :: particles_gat(:)
  real(rp),    pointer     :: particles_send(:)
  type(latyp)              :: particle
  logical(lg), pointer     :: variables_pts(:)
  integer(ip), pointer     :: perm_postprocess(:)
  integer(ip), pointer     :: perm_deposition(:)

  nullify(num_particles_gat) 
  nullify(particles_gat)
  nullify(particles_send)
  nullify(variables_pts)
  nullify(perm_postprocess)
  nullify(perm_deposition)
  
  !----------------------------------------------------------------------
  !
  ! Mesh dependent postprocess
  !
  !----------------------------------------------------------------------

  call output_postprocess_variables(pts_outvar)

  !----------------------------------------------------------------------
  !
  ! Witness
  !
  !----------------------------------------------------------------------

  if( ittyp == ITASK_ENDTIM ) then
     call pts_outwig()
  end if

  !----------------------------------------------------------------------
  !
  ! Particles output
  !
  !----------------------------------------------------------------------
  !
  ! Check what kind of postprocess
  ! KFL_ISPOS = 1: All particles info is required
  ! KFL_IFDEP = 1: Deposited particles info is required
  ! KFL_IFWAL = 1: Particles in moving wall info is required
  ! Initial solution if this is a restart is not written
  !
  kfl_ifpos = 0
  kfl_ifdep = 0
  if( mod(ittim, kfl_oufre_pts) == 0 )                   kfl_ifpos = 1
  if( kfl_injec == 1 )                                   kfl_ifpos = 1
  if( nlagr_going_out_pts > 0 .and. kfl_oudep_pts /= 0 )then
     kfl_ifdep = 1
  end if
  if( ittyp == ITASK_INITIA .and. kfl_rstar == 2 ) then
     kfl_ifpos = 0
     kfl_ifdep = 0
  end if
  !
  ! What type of particle to output
  !
  imaxi =  100000 
  imini = -100000
  if(      kfl_injec == 1 ) then
     imaxi = PTS_PARTICLE_JUST_INJ      ! Only just deposited particles information is gathered
  !else if( kfl_injec == 0 ) then
  !   if( kfl_ifpos == 0 .and. kfl_ifdep == 1 ) then
  !      !   imaxi = PTS_PARTICLE_HITS_WALL  ! Only deposited and vanishing particles information is gathered
  !      !else if ( kfl_ifpos == 0 .and. kfl_lost_moving_pts > 0 ) then
  !      imaxi = PTS_PARTICLE_MOVING_MESH
  !   end if
  end if
  imaxi = min(imaxi,kfl_max_out_exist_state_pts)
  imini = max(imini,kfl_min_out_exist_state_pts)
  if( imaxi < imini )       kfl_ifpos = 0       ! Avoid counting if non of the existence states are requested               
  if( ittim == ittim_last ) kfl_ifpos = 0       ! Last time was already postprocessed
  
  !----------------------------------------------------------------------
  !
  ! MPI_GATHER particles info. Information for postprocess and deposition
  ! are merged in order to avoid extra communication
  !
  !----------------------------------------------------------------------

  if( kfl_ifpos == 1 .or. kfl_ifdep == 1 ) then
     call ker_timeline('INI_OUTPUT')
     !
     ! Save last time we got here
     !
     if( kfl_injec == 1 ) then
        ittim_last = -1
     else
        ittim_last = ittim
     end if
     !
     ! Current time
     !
     if( kfl_injec == 0 ) then
        cutim_pts = cutim
     else
        cutim_pts = cutim-dtime
     end if
     !
     ! Gather
     !
     if( IPARALL ) then
        !
        ! Total number of variables to gather
        !
        call memory_alloca(mem_modul(1:2,modul),'PERM_POSTPROCESS','pts_inivar',perm_postprocess,nvarp_pts)
        call memory_alloca(mem_modul(1:2,modul),'PERM_DEPOSITION', 'pts_inivar',perm_deposition, nvard_pts)
        call memory_alloca(mem_modul(1:2,modul),'VARIABLES_PTS'   ,'pts_inivar',variables_pts,   mvarp_pts)

        nvart_pts = 0
        ivarp     = 0
        ivard     = 0
        do ivari = 1,mvarp_pts
           if(  ( postprocess_var_pts(ivari) .and. kfl_ifpos == 1 ) .or. &
                ( deposition_var_pts(ivari)  .and. kfl_ifdep == 1 ) ) then
              nvart_pts = nvart_pts + 1
          end if
           if( postprocess_var_pts(ivari) .and. kfl_ifpos == 1 ) then
              ivarp                   = ivarp + 1
              perm_postprocess(ivarp) = nvart_pts
              variables_pts(ivari)    = .true.
           end if
           if( deposition_var_pts(ivari)  .and. kfl_ifdep == 1 ) then
              ivard                   = ivard + 1
              perm_deposition(ivard)  = nvart_pts
              variables_pts(ivari)    = .true.
           end if
        end do
        nvarp_pts4 = int(nvarp_pts,4)
        nvard_pts4 = int(nvard_pts,4)
        nvart_pts4 = int(nvart_pts,4)
        num_particles = 0_4

        if( IMASTER ) then
           call memory_alloca(mem_modul(1:2,modul),'NUM_PARTICLES_GAT','pts_output',num_particles_gat,int(npart+1_ip,4_ip),LBOUN=0_4)
        else if( ISLAVE ) then
           num_particles = count( lagrtyp(1:mlagr) % kfl_exist <= imaxi .and. lagrtyp(1:mlagr) % kfl_exist >= imini ) * nvart_pts4
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send,int(num_particles,ip))
           ipars = 0
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist <= imaxi .and. lagrtyp(ilagr) % kfl_exist >= imini ) then
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,variables_pts)
              end if
           end do
        end if
        call PAR_GATHER(num_particles,num_particles_gat)
        if( IMASTER ) then
           num_particles = sum(num_particles_gat)
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_GAT','pts_output',particles_gat,num_particles)
        end if
        call PAR_GATHERV(particles_send,particles_gat,num_particles_gat)
     end if

     !----------------------------------------------------------------------
     !
     ! Postprocess
     !
     !----------------------------------------------------------------------
     
     if( kfl_ifpos == 1 )  then
        if( IMASTER ) then

           ipars = 0
           do ipart = 1,npart
              do ilagr = 1,int(num_particles_gat(ipart),ip) / nvart_pts
                 call pts_parallelization_unpack(particle,ipars,particles_gat,variables_pts,'EXIST')
                 if( particle % kfl_exist <= imaxi .and. particle % kfl_exist >= imini ) then

                     call pts_result_io % write(     particles_gat(perm_postprocess(1)+ipars),&
                                            int(particles_gat(perm_postprocess(2)+ipars),ip),&
                                            int(particles_gat(perm_postprocess(3)+ipars),ip),&                        
                                            int(particles_gat(perm_postprocess(4)+ipars),ip),&                        
                                            particles_gat(perm_postprocess(5:nvarp_pts)+ipars) )
                 end if 
                 ipars = ipars + nvart_pts
              end do
           end do

           call pts_result_io % flush()

        else if( ISEQUEN ) then
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send,nvarp_pts)
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist <= imaxi .and. lagrtyp(ilagr) % kfl_exist >= imini ) then
                 ipars = 0
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,postprocess_var_pts)
                 call pts_result_io % write(     particles_send(1), & 
                                        int(particles_send(2),ip),&
                                        int(particles_send(3),ip),&
                                        int(particles_send(4),ip),&
                                        particles_send(5:nvarp_pts) )
              end if 
           end do
           call memory_deallo(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send)           

           call pts_result_io % flush()

           
        end if

     end if
     
     !----------------------------------------------------------------------
     !
     ! Postprocess deposition
     !
     !----------------------------------------------------------------------
     
     if( kfl_ifdep == 1 ) then

        if( IMASTER ) then

           ipars = 0
           do ipart = 1,npart
              do ilagr = 1,int(num_particles_gat(ipart),ip) / nvart_pts
                 call pts_parallelization_unpack(particle,ipars,particles_gat,variables_pts,'EXIST')
                 if(    particle % kfl_exist == PTS_PARTICLE_HITS_WALL   .or. &
                      & particle % kfl_exist == PTS_PARTICLE_OUTFLOW     .or. &
                      & particle % kfl_exist == PTS_PARTICLE_MOVING_MESH ) then         
                    write(lun_oudep_pts,101) (particles_gat(perm_deposition(ivard)+ipars),ivard=1,nvard_pts)
                 end if
                 ipars = ipars + nvart_pts
              end do 
           end do
           call iofile_flush_unit(lun_oudep_pts)

        else if( ISEQUEN ) then
           
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send,nvard_pts)
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist <= imaxi .and. lagrtyp(ilagr) % kfl_exist >= imini ) then
                 if(    lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_HITS_WALL   .or. &
                      & lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_OUTFLOW     .or. &
                      & lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_MOVING_MESH ) then
                    ipars = 0
                    call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,deposition_var_pts)
                    write(lun_oudep_pts,101) (particles_send(ivard),ivard=1,nvard_pts)
                 end if
              end if
           end do
           call memory_deallo(mem_modul(1:2,modul),'PARTICLES_SEND','pts_output',particles_send)
           call iofile_flush_unit(lun_oudep_pts)

        end if

     end if

     call ker_timeline('END_OUTPUT')

  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess surface deposition 
  !
  !----------------------------------------------------------------------

  if( kfl_ifdep == 1 .and. kfl_depos_surface_pts == 1 .and. nlagr_hits_wall_pts > 0 ) then

     allocate( depos_surface(ntyla_pts+1) )
     call pts_deposition_surface(depos_surface)

     if( INOTSLAVE ) then
        write(lun_depsu_pts,101) cutim_pts,depos_surface(ntyla_pts+1),depos_surface(1:ntyla_pts)
        call iofile_flush_unit(lun_depsu_pts)
     end if

     deallocate( depos_surface )

  end if
  !
  ! Deallocate memory
  !
  call memory_deallo(mem_modul(1:2,modul),'NUM_PARTICLES_GAT' ,'pts_output',num_particles_gat)
  call memory_deallo(mem_modul(1:2,modul),'PARTICLES_SEND'    ,'pts_output',particles_send)
  call memory_deallo(mem_modul(1:2,modul),'PARTICLES_GAT'     ,'pts_output',particles_gat)
  call memory_deallo(mem_modul(1:2,modul),'PERM_POSTPROCESS'  ,'pts_output',perm_postprocess)
  call memory_deallo(mem_modul(1:2,modul),'PERM_DEPOSITION'   ,'pts_output',perm_deposition)
  call memory_deallo(mem_modul(1:2,modul),'VARIABLES_PTS'     ,'pts_output',variables_pts)
  !
  ! Formats
  !
101 format(es16.8e3,400(',',es16.8e3))

end subroutine pts_output


subroutine pts_deposition_header()

  use def_master
  use def_partis
  implicit none

  if( nvard_pts > 0 ) then
     write(lun_oudep_pts,1) postprocess_name_pts(deposition_list_pts(1:nvard_pts))
  end if
  
1 format(100(a5,','))
  
end subroutine pts_deposition_header
