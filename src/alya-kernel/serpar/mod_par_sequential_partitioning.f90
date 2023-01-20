!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_sequential_partitioning.f90
!> @author  houzeaux
!> @date    2018-10-09
!> @brief   Sequential partitioning
!> @details Tools for sequential partitioning
!-----------------------------------------------------------------------

module mod_par_sequential_partitioning

  use def_parame
  use def_master
  use def_domain 
  use def_parall
  use mod_memory,                    only : memory_alloca
  use mod_reabcs,                    only : reabcs_seq
  use mod_partition_sfc,             only : partition_sfc
  use mod_redistribute,              only : par_redistribute
  use mod_redistribute,              only : gather_to_master 
  use mod_parall,                    only : PAR_METIS4
  use mod_parall,                    only : PAR_SFC
  use mod_parall,                    only : PAR_ORIENTED_BIN
  use mod_parall,                    only : PAR_PARALLEL_PARTITION
  use mod_communications,            only : PAR_BARRIER
  use mod_communications,            only : PAR_BROADCAST
  use def_master,                    only : kfl_paral  
  use def_domain,                    only : ndime
  use def_master,                    only : new_periodicity
  use mod_parall,                    only : PAR_COMM_MY_CODE_WM
  use mod_parall,                    only : PAR_WORLD_SIZE
  use mod_parall,                    only : par_memor
  use def_parall,                    only : boxes_coarse_par 
  use def_parall,                    only : boxes_fine_par 
  use mod_par_interface_exchange,    only : par_interface_exchange
  use mod_parall,                    only : PAR_COMM_MY_CODE_ARRAY
  use mod_par_additional_arrays,     only : par_ordered_exchange_update
  use mod_domain,                    only : domain_memory_allocate     
  use mod_domain,                    only : domain_memory_deallocate     

  implicit none
  
  public :: par_sequential_partitioning
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-09
  !> @brief   Main subroutine
  !> @details Perform sequential partitioning
  !> 
  !-----------------------------------------------------------------------

  subroutine par_sequential_partitioning()

  integer(ip)          :: ichec,ipart,iauxi,dummi,npart_sfc
  integer(ip)          :: npoin_tmp, nelem_tmp, nboun_tmp
  real(rp)             :: time1,time2,time3,time4,time5
  character(100)       :: messa_integ
  real(rp),    pointer :: lenti_(:,:) 
  integer(ip), pointer :: lweig_(:)   
  integer(4),  pointer :: nelem_part_master4_(:) 
  integer(4),  pointer :: npoin_part_master4_(:) 

  nullify(lenti_,lweig_,nelem_part_master4_,npoin_part_master4_)

  call PAR_BARRIER()
  
  !-------------------------------------------------------------------
  !
  ! Compute graphs, needed for sequential partitioning and also groups
  ! if the sequential frontal method has been selected
  !
  !-------------------------------------------------------------------

  if( IMASTER .and. .not. READ_AND_RUN() ) then
     !
     ! Compute graph
     !
     call domgra(1_ip)                          ! PELPO and LELPO, R_DOM, C_DOM
     !
     ! Groups for deflated (only in sequential mode)
     !
     call grodom(1_ip)
  end if

  if( INOTMASTER ) then
       
       !-------------------------------------------------------------------
       !
       ! SFC
       !
       !-------------------------------------------------------------------

       if( kfl_partition_par == PAR_SFC ) then
          npart_sfc = PAR_WORLD_SIZE - 1_ip
          call runend('SFRC PARTITIONING SHOULD BE USED WITH PARALLEL EXECUTION MODE')
          !call par_redistribute(npart_sfc,lenti_,lweig_,nenti_,nelem_part_master4_,npoin_part_master4_)
          !call partition_sfc(lepar_par,npart_sfc,lenti_,nenti_,ndime,lweig_,boxes_coarse_par,boxes_fine_par,PAR_COMM_=PAR_COMM_MY_CODE_WM)
          !call gather_to_master(lepar_par,npart)
       end if

    end if

    call cputim(time1)
    
    if( IMASTER .and. .not. READ_AND_RUN() ) then

       !-------------------------------------------------------------------
       !
       ! Master
       !
       !-------------------------------------------------------------------

       call cputim(time1)
       call par_livinf(1_ip,' ',dummi)
       !
       ! Prepare arrays needed for efficient reinitialization in cut nodes,  LEFPO & PEFPO in master
       !
       call cputim(time2)
       !
       ! Partition mesh
       !
       call cputim(time4)

       npart_sfc = PAR_WORLD_SIZE - 1_ip
       call par_partit()
       
       call cputim(time5)
       !
       ! Compute arrays
       !
       call par_arrays()
       !
       ! Send data to slaves
       !
       call par_livinf(6_ip,' ',dummi)

       call cputim(time4)
       do kfl_desti_par = 1,npart_par

          iauxi = npart_par / 10_ip
          messa_integ='    S/R  '//trim(intost(kfl_desti_par))//' - TOT '//trim(intost(npart_par))

          call par_sendat(two)            ! Send/receive data from readim and reastr
          call par_sendat(eight)          ! Send/receive data from readim and reastr
          call par_sendat(five)           ! Send/receive data from partit

          if (iauxi > 0) then
             if (modulo(kfl_desti_par,iauxi)==0) then
                call par_livinf(20_ip,messa_integ,dummi)
             end if
          end if

       end do
       
       if(associated(lenti_))              deallocate(lenti_)
       if(associated(lweig_))              deallocate(lweig_)
       if(associated(nelem_part_master4_)) deallocate(nelem_part_master4_)
       if(associated(npoin_part_master4_)) deallocate(npoin_part_master4_)

       messa_integ='  SEND/REC  '//trim(intost(npart_par))//' - TOTAL '//trim(intost(npart_par))
       call par_livinf(20_ip,messa_integ,dummi)

       call par_livinf(20_ip,'     COORD, LTYPE, LNODS, LTYPB, LNODB... ',dummi)
       call par_sengeo(1_ip)              
       call par_livinf(20_ip,'     LELBO... ',dummi)
       call par_sengeo(2_ip)       
       call par_livinf(20_ip,'     KFL_FIELD, XFIEL, TIME_FIELD,... ',dummi)
       call par_sengeo(3_ip)    
       call par_livinf(20_ip,'     SET DATA... ',dummi)
       call par_senset()         
       call par_livinf(20_ip,'     BC DATA... ',dummi)
       call par_senbcs()
       call par_livinf(20_ip,'     COMMUNICATION ARRAYS... ',dummi)
       call par_sencom()           
       call par_livinf(20_ip,'  SEND DATA TO SLAVES... DONE.',dummi)

       call cputim(time5)
       !
       ! Barrier
       !
       do ipart = 1,nproc_par-1,nsire_par
          call PAR_BARRIER()
       end do
       !
       ! Local ordering of nodes
       !
       npoi1 =  0
       npoi2 =  0
       npoi3 = -1
       npoin =  0
       nelem =  0
       nboun =  0
       !
       ! Writes in partition file
       !
       call par_livinf(1000_ip,' ',dummi)
       if( PART_AND_WRITE() ) call par_outprt()
       call cputim(time3)
       cpu_paral(4) = time3 - time2

    else if( ISLAVE .and. ( PART_AND_RUN() .or. READ_AND_RUN() ) ) then

       !-------------------------------------------------------------------
       !
       ! Slaves
       !
       !-------------------------------------------------------------------
       !
       ! Receive (kfl_ptask==1) of read (kfl_ptask==2) data with barrier
       !
       do ipart = 1,nproc_par-1,nsire_par
          if( kfl_paral >= ipart .and. kfl_paral < ipart+nsire_par ) then
             kfl_desti_par = 0
             call par_sendat(two)                    ! Send/receive data from readim and reastr
             call par_sendat(eight)                  ! Send/receive data from readim and reastr
             call par_sendat(five)                   ! Send/receive data from partit
             !
             ! Local ordering of nodes
             !
             npoi1 = lni                             ! Interior nodes: ipoin=1,npoi1
             npoi2 = slfbo                           ! begin of boundary nodes: ipoin = npoi2,npoi3
             npoi3 = lnb + slfbo - 1                 ! End of own boundary nodes
             !
             ! Continue sending
             !
             call domain_memory_allocate('GEOMETRY') ! Allocate memory for geometrical arrays     
             call par_sengeo(1_ip)                   ! Receive geometry: COORD, LTYPE, LNODS, LTYPB, LNODB, IB...
             call par_sengeo(2_ip)                   ! Receive geometry: LELBO
             call par_sengeo(3_ip)                   ! Receive geometry: KFL_FIELD, XFIEL, TIME_FIELD
             call par_senset()                       ! Receive set data
             call par_senbcs()                       ! Receive boundary conditions
             call par_sencom()                       ! Receive communication arrays
          end if
          call PAR_BARRIER()

       end do
       !
       ! Checkpoint
       !
       if( READ_AND_RUN() ) call par_chkpoi(ichec)

    else if( IMASTER .and. READ_AND_RUN() ) then

       !-------------------------------------------------------------------
       !
       ! Master in read-and-run mode
       !
       !-------------------------------------------------------------------

       if (kfl_modul(id_levels)/=0) call runend('par_prepro: read_and_run not ready with level set')
       !
       ! Open minimum graph to avoid null pointers when calling solvers
       !
       call memory_alloca(par_memor,'R_DOM','par_memory',r_dom,1_ip)
       call memory_alloca(par_memor,'C_DOM','par_memory',c_dom,1_ip)
       !
       ! Local ordering of nodes
       !
       npoi1 =  0
       npoi2 =  0
       npoi3 = -1
       npoin =  0
       nelem =  0
       nboun =  0 ! Some time ago it didn't work, let's see now
       !
       ! Barrier
       !
       do ipart = 1,nproc_par-1,nsire_par
          call PAR_BARRIER()
       end do
       !
       ! Communicators
       !
       call par_sencom()
       !
       ! Master reads from partition file
       !
       call par_livinf(6_ip,' ',dummi)
       call par_outprt()
       !
       ! Checkpoint
       !
       call par_chkpoi(ichec)
       call par_livinf(10_ip,' ',ichec)

    end if
    
  end subroutine par_sequential_partitioning
  
end module mod_par_sequential_partitioning
!> @}
