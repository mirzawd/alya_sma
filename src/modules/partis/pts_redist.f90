!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_redist.f90
!> @author  houzeaux
!> @date    2020-02-04
!> @brief   Redistribution of particles
!> @details redistribute particles
!> @} 
!-----------------------------------------------------------------------

subroutine pts_redist()

  use def_master
  use def_domain
  use def_partis
  use mod_redistribution, only : commd_nelem
  use mod_redistribution
  use mod_pts_parallelization
  use mod_communications
  use mod_memory
  use mod_messages
  use mod_pts_arrays, only : pts_arrays
  implicit none
  integer(ip)            :: ilagr_local,ilagr,ielem,ipars
  integer(ip)            :: iparr,icror,ifoun,ineig,ii,nlagr_last
  integer(ip)            :: new_size,kelem,nelem_max,kk,nn
  integer(ip),  pointer  :: counte(:) 
  integer(ip),  pointer  :: perme(:) 
  real(rp),     pointer  :: send_pts(:)
  real(rp),     pointer  :: recv_pts(:)
  integer(ip),  pointer  :: nlagr_send(:)
  integer(ip),  pointer  :: nlagr_recv(:)
  type(i1p),    pointer  :: list_particles(:)

  nullify(send_pts)
  nullify(recv_pts)
  nullify(nlagr_send)
  nullify(nlagr_recv)
  nullify(counte)
  nullify(perme)
  nullify(list_particles)
  particles_sent = 0
  particles_recv = 0

  if( INOTMASTER ) then
     !
     ! NELEM is the number of element of the new mesh...
     ! Try to guess the old number of elements
     !
     nelem_max = 0
     do ii = 1,commd_nelem % lsend_dim
        ielem = commd_nelem % lsend_perm(ii)
        nelem_max = max(nelem_max,ielem)
     end do
     !
     ! Count number of particles per elements
     !     
     call memory_alloca(mem_modul(1:2,modul),'COUNTE'        ,'pts_solite',counte,nelem_max) ! should be old nelem
     call memory_alloca(mem_modul(1:2,modul),'LIST_PARTICLES','pts_solite',list_particles,nelem_max) 
     do ilagr_local = 1,nlagr_local_pts 
        ilagr         = permu_nlagr_pts(ilagr_local)
        ielem         = int(lagrtyp(ilagr) % ielem,ip)
        counte(ielem) = counte(ielem)+1
     end do
     !
     ! Store particles in elements
     !
     do ielem = 1,nelem_max
        call memory_alloca(mem_modul(1:2,modul),'LIST_PARTICLES % L','pts_solite',list_particles(ielem) % l,counte(ielem))
        counte(ielem) = 0
     end do
     do ilagr_local = 1,nlagr_local_pts
        ilagr         = permu_nlagr_pts(ilagr_local)
        ielem         = int(lagrtyp(ilagr) % ielem,ip)
        counte(ielem) = counte(ielem)+1
        list_particles(ielem) % l(counte(ielem)) = ilagr
     end do
     !
     ! # particles to send to each neighbor
     !
     allocate(nlagr_send(commd_nelem % nneig))
     allocate(nlagr_recv(commd_nelem % nneig))

     do ineig = 1,commd_nelem % nneig
        nlagr_send(ineig) = 0
        do ii = commd_nelem % lsend_size(ineig),commd_nelem % lsend_size(ineig+1)-1
           ielem = commd_nelem % lsend_perm(ii)
           nlagr_send(ineig) = nlagr_send(ineig) + counte(ielem)
        end do
        nlagr_send(ineig) = nlagr_send(ineig)
     end do
     !
     ! # particles to receive from each neighbor
     !
     do ineig = 1,commd_nelem % nneig
        call PAR_SEND_RECEIVE(1_ip,1_ip,nlagr_send(ineig:ineig),nlagr_recv(ineig:ineig),'IN MY ZONE',commd_nelem % neights(ineig) )
     end do
     !
     ! Possibly reallocate
     !
     do ineig = 1,commd_nelem % nneig
        if( commd_nelem % neights(ineig) /= kfl_paral ) then
           particles_sent = particles_sent + nlagr_send(ineig) 
           particles_recv = particles_recv + nlagr_recv(ineig)
        end if
     end do
  end if
  
  nn = particles_sent 
  call PAR_SUM(nn)
  call messages_live('PARTIS: MIGRATING '//trim(intost(nn))//' PARTICLES')
  
  if( INOTMASTER ) then
     
     if( particles_recv > nlagr_free_pts ) then
        new_size = int(1.2_rp*real(mlagr+particles_recv-nlagr_free_pts,rp),ip)
        call pts_reallocate(new_size)
        if( size(permu_nlagr_pts) < mlagr ) then
           call memory_resize(mem_modul(1:2,modul),'PERMU_NLAGR_PTS','pts_solite',permu_nlagr_pts,mlagr)
        end if
     end if
     !
     ! Total size
     !
     do ineig = 1,commd_nelem % nneig
        if( commd_nelem % neights(ineig) /= kfl_paral ) then
           nlagr_send(ineig) = nlagr_send(ineig) * number_migrated_variables_pts
           nlagr_recv(ineig) = nlagr_recv(ineig) * number_migrated_variables_pts
        end if
     end do
     !
     ! Redistribute particles
     !
     nlagr_last = 0

     do ineig = 1,commd_nelem % nneig

        if( commd_nelem % neights(ineig) /= kfl_paral ) then

           call memory_alloca(mem_modul(1:2,modul),'SEND_PTS','pts_solite',send_pts,nlagr_send(ineig),'DO_NOT_INITIALIZE')
           call memory_alloca(mem_modul(1:2,modul),'RECV_PTS','pts_solite',recv_pts,nlagr_recv(ineig),'DO_NOT_INITIALIZE')
           !
           ! Pack particles
           !
           ipars = 0
           do ii = commd_nelem % lsend_size(ineig),commd_nelem % lsend_size(ineig+1)-1
              ielem = commd_nelem % lsend_perm(ii)
              do ilagr_local = 1,counte(ielem)
                 ilagr                  = list_particles(ielem) % l(ilagr_local)
                 lagrtyp(ilagr) % ielem = ii+1-commd_nelem % lsend_size(ineig)
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,send_pts,migrated_variables_pts)             
                 lagrtyp(ilagr) % kfl_exist   = 0          ! Take off particle from my list
              end do
           end do
           !
           ! Send particles
           !
           call PAR_SEND_RECEIVE(send_pts,recv_pts,'IN MY ZONE',commd_nelem % neights(ineig) )
           ! 
           ! Allocate new particles
           !
           ilagr      = 1
           iparr      = 0
           !
           ! Unpack particles and look for a free space in LAGRTYP to save particle
           ! If there is no space, allocate more memory
           ! 
           do icror = 1,nlagr_recv(ineig)/number_migrated_variables_pts
 
              ifoun = 0
              ilagr = nlagr_last
              loop_find_position: do while( ilagr < mlagr )
                 ilagr = ilagr + 1
                 if( lagrtyp(ilagr) % kfl_exist == 0 ) then
                    ifoun = 1
                    exit loop_find_position
                 end if
              end do loop_find_position
              nlagr_last = ilagr
              if( ifoun == 0 ) then
                 print*,'we are in big troubles'
                 call runend('PTS_SOLITE: WE ARE IN TROUBLE!')
              end if
              call lagrtyp(ilagr) % init()           
              lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_EXISTS
              call pts_parallelization_unpack(lagrtyp(ilagr),iparr,recv_pts,migrated_variables_pts)
              kk                        = lagrtyp(ilagr) % ielem-1+commd_nelem % lrecv_size(ineig)
              kelem                     = commd_nelem % lrecv_perm(kk)
              lagrtyp(ilagr) % ielem    = kelem
              lagrtyp(ilagr) % mpi_rank = kfl_paral
           end do
           !
           ! Deallocate 
           !
           call memory_deallo(mem_modul(1:2,modul),'SEND_PTS','pts_solite',send_pts)
           call memory_deallo(mem_modul(1:2,modul),'RECV_PTS','pts_solite',recv_pts)

        else

           do ii = commd_nelem % lsend_size(ineig),commd_nelem % lsend_size(ineig+1)-1
              ielem = commd_nelem % lsend_perm(ii)
              do ilagr_local = 1,counte(ielem)
                 ilagr                  = list_particles(ielem) % l(ilagr_local)
                 lagrtyp(ilagr) % ielem = ii+1-commd_nelem % lsend_size(ineig)                    
                 kk                     = lagrtyp(ilagr) % ielem-1+commd_nelem % lrecv_size(ineig)
                 kelem                  = commd_nelem % lrecv_perm(kk)
                 lagrtyp(ilagr) % ielem = kelem
                 !print*,'poposssssssssssssssssssssop=',ilagr
              end do
           end do

        end if
        
     end do

     call memory_alloca(mem_modul(1:2,modul),'COUNTE'        ,'pts_solite',counte,nelem) 
     call memory_alloca(mem_modul(1:2,modul),'LIST_PARTICLES','pts_solite',list_particles,commd_nelem % nneig)
     !
     ! Recompute permutation
     !
     call pts_compute_permutation() 

  end if

  call pts_arrays('REDISTRIBUTE')
  call redistribution_array(kfl_fixbo_pts           ,'NBOUN',MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_PTS')
  call redistribution_array(bvnat_pts               ,'NBOUN',MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVNAT_PTS')
  call redistribution_array(kfl_fixno_walld_slip_pts,'NPOIN',MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_WALLD_SLIP_PTS')

end subroutine pts_redist
