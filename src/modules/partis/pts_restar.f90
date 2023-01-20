!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @name    Partis restart
!> @file    pts_restar.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   This routine writes or reads restart file
!> @details Restart file. Write and read the same variables as for
!>          the migration of particles through subdomains.
!> @} 
!------------------------------------------------------------------------

subroutine pts_restar(itask)
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_partis
  use mod_memory
  use mod_postpr
  use mod_communications,      only : PAR_GATHER
  use mod_communications,      only : PAR_GATHERV
  use mod_communications,      only : PAR_SUM
  use mod_communications,      only : PAR_SEND_RECEIVE
  use mod_communications,      only : PAR_SEND
  use mod_communications,      only : PAR_RECEIVE
  use mod_pts_parallelization, only : pts_parallelization_pack
  use mod_pts_parallelization, only : pts_parallelization_unpack
  use mod_messages,            only : livinf
  use mod_pts_arrays,          only : pts_arrays
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipart,ilagr,ipars,isize,iinj
  integer(ip)             :: nlagr_0,nvar2,kfl_reawr_tmp,ivari
  integer(4)              :: num_particles4,nvar2_4
  real(rp)                :: dummr_recv(2)
  real(rp)                :: dummr_send(2)
  integer(4),  pointer    :: num_particles_gat4(:)
  real(rp),    pointer    :: particles_gat(:)
  real(rp),    pointer    :: particles_send(:)
  !
  ! Check if this is a preliminary or restart
  !
  !all respre(itask,kfl_gores)
  !if( kfl_gores == 0 ) return

  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then
     call pts_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call pts_arrays('WRITE RESTART')
  end if
  
  nullify(num_particles_gat4)
  nullify(particles_gat)
  nullify(particles_send)
  nvar2   = number_migrated_variables_pts
  nvar2_4 = int(nvar2,4)
  
  
  select case ( itask )

  case ( ITASK_READ_RESTART )

     !-------------------------------------------------------------------
     !
     ! Read continue restart file
     !
     !-------------------------------------------------------------------

     !
     ! Open file
     !
     kfl_reawr_tmp = kfl_reawr
     call PAR_SUM(kfl_reawr_tmp)
     !
     ! File does not exist: return in normal mode
     !
     if( kfl_reawr_tmp < 0 ) return

     if( IMASTER ) then
        !
        ! Read restart file
        ! IPARS .............................. Total number of particles * NVAR2
        ! num_particles_gat4(IPART) / NVAR2 ... Number of particles to read for subdomain IPART
        !           
        call memory_alloca(mem_modul(1:2,modul),'NUM_PARTICLES_GAT4','pts_restar',num_particles_gat4,npart,'DO_NOT_INITIALIZE')

        read(momod(modul) % lun_rstar) nlacc_pts
        do iinj= 1, pts_minj
           read(momod(modul) % lun_rstar) injection_pts(iinj) % time_cumulative
        enddo
        read(momod(modul) % lun_rstar) (num_particles_gat4(ipart),ipart=1,npart)

        ipars = sum(num_particles_gat4)
        call livinf(-9_ip,'READ LAGRANGIAN PARTICLES= ',ipars / nvar2) 
        call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_restar',particles_send,ipars,'DO_NOT_INITIALIZE')

        ipars = 0
        do ipart = 1,npart
           do ilagr = 1,num_particles_gat4(ipart) / nvar2
              read(momod(modul) % lun_rstar) (particles_send(ipars+ivari),ivari=1,nvar2)
              ipars = ipars + nvar2
           end do
        end do

        ipars = 1
        do ipart = 1,npart              
           isize = num_particles_gat4(ipart)              
           call PAR_SEND(nlacc_pts,'IN MY CODE',ipart)
           do iinj= 1, pts_minj
              call PAR_SEND(injection_pts(iinj) % time_cumulative,'IN MY CODE',ipart)
           enddo
           call PAR_SEND(isize    ,'IN MY CODE',ipart)
           if( isize > 0 ) then
              call PAR_SEND_RECEIVE(isize,0_ip,particles_send(ipars:),dummr_recv,'IN MY CODE',ipart)
           end if              
           ipars = ipars + isize
        end do
        
     else if( ISEQUEN ) then
        !
        ! Sequen read file
        !        
        read(momod(modul) % lun_rstar) nlacc_pts
        do iinj= 1, pts_minj
           read(momod(modul) % lun_rstar) injection_pts(iinj) % time_cumulative
        enddo
        read(momod(modul) % lun_rstar) isize
        call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_restar',particles_send,nvar2)

        do ilagr = 1,isize
           lagrtyp(ilagr) % kfl_exist = -1
           ipars = 0
           read(momod(modul) % lun_rstar) (particles_send(ivari),ivari=1,nvar2)
           call pts_parallelization_unpack(lagrtyp(ilagr),ipars,particles_send,migrated_variables_pts)
           if (kfl_rstar == 1) then
               !
               ! Deal with initialization 
               !
               lagrtyp(ilagr) % t        = 0.0_rp
               lagrtyp(ilagr) % t_inject = 0.0_rp
           endif
        end do

     else if( ISLAVE ) then
        !
        ! Slave receive particle info from master
        !
        ipart = 0
        call PAR_RECEIVE(nlacc_pts,'IN MY CODE',ipart)
        do iinj= 1, pts_minj
           call PAR_RECEIVE(injection_pts(iinj) % time_cumulative,'IN MY CODE',ipart)
        enddo
        call PAR_RECEIVE(isize    ,'IN MY CODE',ipart)

        nlagr_0 = isize / nvar2
        
        if( nlagr_0 > mlagr ) call pts_reallocate(nlagr_0)

        if( isize > 0 ) then
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_restar',particles_send,isize,'DO_NOT_INITIALIZE')
           call PAR_SEND_RECEIVE(0_ip,isize,dummr_send,particles_send,'IN MY CODE',ipart)
           ipars = 0
           do ilagr = 1,nlagr_0
              lagrtyp(ilagr) % kfl_exist = -1
              call pts_parallelization_unpack(lagrtyp(ilagr),ipars,particles_send,migrated_variables_pts)
              if (kfl_rstar == 1) then
                  !
                  ! Deal with initialization 
                  !
                  lagrtyp(ilagr) % t        = 0.0_rp
                  lagrtyp(ilagr) % t_inject = 0.0_rp
              endif
           end do              
        end if 

     end if

  case ( ITASK_WRITE_RESTART )

     !-------------------------------------------------------------------
     !
     ! Write restart file
     !
     !-------------------------------------------------------------------

     if( IPARALL ) then
        num_particles4 = 0_4
        if( IMASTER ) then
           call memory_alloca(mem_modul(1:2,modul),'NUM_PARTICLES_GAT4','pts_restar',num_particles_gat4,int(npart+1_ip,4),LBOUN=0_4)
        else if( ISLAVE ) then
           num_particles4 = count( lagrtyp(1:mlagr) % kfl_exist == -1 ) * nvar2_4
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_restar',particles_send,int(num_particles4,ip))
           ipars = 0
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist == -1 ) then
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,migrated_variables_pts)
              end if 
           end do
        end if
        call PAR_GATHER(num_particles4,num_particles_gat4)
        if( IMASTER ) then
           num_particles4 = sum(num_particles_gat4)
           call memory_alloca(mem_modul(1:2,modul),'PARTICLES_GAT','pts_restar',particles_gat,num_particles4)
        end if
        !
        ! MPI_Gatherv  
        !
        call PAR_GATHERV(particles_send,particles_gat,num_particles_gat4)
     end if
     
     if( IMASTER ) then
        !
        ! Only write existing particles with kfl_exist = -1
        ! PARIG(IPART) ......... Length of data received by slave IPART
        ! PARIG(IPART)/NVAR2 ... Number of existing particles received by slave IPART
        !
        write(momod(modul) % lun_rstar) nlacc_pts
        do iinj= 1, pts_minj
           write(momod(modul) % lun_rstar) injection_pts(iinj) % time_cumulative
        enddo
        write(momod(modul) % lun_rstar) (num_particles_gat4(ipart),ipart=1,npart)
        ipars = 0
        do ipart = 1,npart
           do ilagr = 1,num_particles_gat4(ipart) / nvar2
              write(momod(modul) % lun_rstar) (particles_gat(ipars+ivari),ivari=1,nvar2)
              ipars = ipars + nvar2
           end do
        end do

     else if( ISEQUEN ) then

        call memory_alloca(mem_modul(1:2,modul),'PARTICLES_SEND','pts_restar',particles_send,nvar2)
        isize = 0
        do ilagr = 1,mlagr
           if( lagrtyp(ilagr) % kfl_exist == -1 ) isize = isize + 1
        end do
        write(momod(modul) % lun_rstar) nlacc_pts
        do iinj= 1, pts_minj
           write(momod(modul) % lun_rstar) injection_pts(iinj) % time_cumulative
        enddo
        write(momod(modul) % lun_rstar) isize
        do ilagr = 1,mlagr
           if( lagrtyp(ilagr) % kfl_exist == -1 ) then
              ipars = 0
              call pts_parallelization_pack(lagrtyp(ilagr),ipars,particles_send,migrated_variables_pts)
              write(momod(modul) % lun_rstar) (particles_send(ivari),ivari=1,nvar2)
           end if
        end do

     end if

  end select

  call memory_deallo(mem_modul(1:2,modul),'NUM_PARTICLES_GAT4','pts_restar',num_particles_gat4)
  call memory_deallo(mem_modul(1:2,modul),'PARTICLES_SEND'    ,'pts_restar',particles_send)
  call memory_deallo(mem_modul(1:2,modul),'PARTICLES_GAT'     ,'pts_restar',particles_gat)

end subroutine pts_restar
