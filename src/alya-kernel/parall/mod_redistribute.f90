!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!





module mod_redistribute

   use mod_communications,  only : PAR_COMM_SPLIT
   use mod_parall,          only : PAR_COMM_SFC
   use mod_parall,          only : PAR_COMM_SFC_WM
   use mod_parall,          only : PAR_MY_SFC_RANK_WM
   use mod_parall,          only : PAR_MY_CODE_RANK
   use mod_parall,          only : PAR_COMM_MPIO, PAR_COMM_MPIO_WM, PAR_COMM_MPIO_RANK_WM, PAR_COMM_MPIO_WM_SIZE
   use mod_parall,          only : par_memor
   use def_parall,          only : boxes_fine_par
   use def_parall,          only : boxes_coarse_par
   use def_parall,          only : sfc_criteria
   use def_parall,          only : sfc_check
   use mod_memory,          only : memory_alloca, memory_deallo
   use def_domain,          only : nnode, ndime, mnode, nelem, npoin, ngaus, nelty
   use def_domain,          only : lnods,ltype,coord
   use def_master,          only : IPARSLAVE, IMASTER,INOTMASTER
   use def_kintyp,          only : ip, rp, lg, i1p, r2p
   use mod_communications,  only : PAR_SUM, PAR_BARRIER
   use mod_communications,  only : PAR_SEND_RECEIVE
   use mod_communications,  only : PAR_START_NON_BLOCKING_COMM
   use mod_communications,  only : PAR_END_NON_BLOCKING_COMM
   use mod_communications,  only : PAR_SET_NON_BLOCKING_COMM_NUMBER
   use mod_communications,  only : PAR_COMM_RANK_AND_SIZE, PAR_MAX
   use mod_mpio_par_configure,  only : dboxc_p

   implicit none


   ! Partition criteria: 0=nodes, 1=elements, 2=weigh.elem.(gauss points)
                                                        !                     3=weigh. elem. (entries)
   !rick: falta comprobar que CRITERIA es 0,1,2 y deixar de comprobar en cada lloc
   private

   integer(ip)             :: dboxc(3)        ! #bin boxes per direction in the coarse bin
   integer(ip)             :: nboxc           ! #boxes of the coarse bin

   integer(ip)             :: nelem_part      ! Number of elements assigned to the partitioner !rick: treure?
   integer(ip)             :: npoin_part      ! Number of nodes assigned to the partitioner !rick: treure?
   integer(ip)             :: npoin_aver      ! Number of points of all partitions except last one !rick: treure
   integer(ip)             :: ipoi0           ! Starting index to go from your local ipoin to the global one !rick: treure
   real(rp),    pointer    :: coord_par(:,:)  ! coordinates of each partitioner slave
   real(rp),    pointer    :: coord_loc(:,:)  ! coordinates of each partitioner slave with local numbering
   integer(ip), pointer    :: ltype_par(:)    ! type of element of each partitioner slave
   integer(ip), pointer    :: lnods_par(:,:)  ! interior element connectivity for each partitioner slave
   integer(ip), pointer    :: lnods_loc(:,:)  ! interior element connectivity for each partitioner slave with local numbering
   integer(ip), pointer    :: invdict(:)      ! Dictionary to go from global to local numbering invdict(old) = new

   integer(4),  pointer    :: nelem_part_master4(:)  ! Number of elements sent to a partition (only the master has it)
   integer(4),  pointer    :: npoin_part_master4(:)  ! Number of points sent to a partition (only the master has it)
   integer(4)              :: comm_rank, comm_size
   real(rp), pointer       :: lenti(:,:)
   integer(ip),pointer     :: lweig(:)
   integer(ip)             :: nenti
   !
   ! Public functions
   !
   public :: par_redistribute
   public :: gather_to_master
   public :: par_visualization_sfc

contains

   !----------------------------------------------------------------------
   !
   !> @author  Juan Cajas & Ricard Borrell
   !> @date    16/03/2018
   !> @brief   Redistribute mesh before SFC
   !> @details
   !
   !----------------------------------------------------------------------
   subroutine par_redistribute(npart,lenti_,lweig_,nenti_,nelem_part_master4_,npoin_part_master4_)

      implicit none
      integer(ip),           intent(inout)    :: npart        ! #partitions
      integer(ip),           intent(inout)    :: nenti_
      real(rp),    pointer,  intent(inout)    :: lenti_(:,:)
      integer(ip), pointer,  intent(inout)    :: lweig_(:)
      integer(4),  pointer,  intent(inout)    :: nelem_part_master4_(:)
      integer(4),  pointer,  intent(inout)    :: npoin_part_master4_(:)

      !
      ! Make arguments vaiables of de module
      !


#ifndef MPI_OFF
      !
      ! Define slaves communicator
      !
      call par_slaves_communicator()
      !
      ! Distribute the mesh among the slaves
      !
      call par_distribute_mesh()
      !
      ! Exchange coordinates needed among slaves
      !
      call par_exchange_coordinates()
      !
      ! Prepare inputs for sfc partition
      !
      call par_sfc_input()
      !
      ! Deallocate module pointers
      !
      call par_deallocate()
#endif
      !
      ! Make arguments vaiables of de module
      !
      nenti_ =  nenti
      lenti_ => lenti
      lweig_ => lweig
      nelem_part_master4_ => nelem_part_master4
      npoin_part_master4_ => npoin_part_master4

    end subroutine par_redistribute

   !
   ! Define slaves communicator
   !
   subroutine par_slaves_communicator()

      use mod_communications, only : PAR_COMM_SPLIT, PAR_BROADCAST
      use mod_parall,         only : PAR_CODE_SIZE

      implicit none
      integer(ip)           :: kfl_paral_proc_node,ndims, auxi
!      integer(ip)           :: ipart
      integer(4_ip)         :: icolor,iproc


#ifndef NDIMEPAR
      if( IMASTER )    ndims = ndime
      call PAR_BROADCAST(ndims,'IN THE WORLD')
      if( INOTMASTER ) ndime = ndims
#endif

      !
      ! Set number of slaves according to hilbert restriction
      !
      if( boxes_fine_par(1) == 0 .and. boxes_coarse_par(1) == 0 ) then
         if( ndime == 3_ip ) then
            auxi       = int(log(real(PAR_CODE_SIZE,rp))/log(8.0_rp),ip)
            nboxc      = 8_ip**auxi
            dboxc(1:3) = 2_ip**auxi
         else if( ndime == 2_ip ) then
            auxi       = int(log(real(PAR_CODE_SIZE,rp))/log(4.0_rp),ip)
            nboxc      = 4_ip**auxi
            dboxc(1:3) = 2**auxi
         else
            call runend("par_slaves_communicator: ndime/=2 and ndime/=3")
         endif
      else
         nboxc = boxes_coarse_par(1)**ndime
      endif
      !
      ! Create the slaves communicator
      !
      kfl_paral_proc_node = 1_ip

      IPARSLAVE=.TRUE.
      if(IMASTER) IPARSLAVE=.FALSE.

      !IPARSLAVE = .FALSE.
      !do ipart=0,nboxc - 1_ip
      !   if (PAR_MY_CODE_RANK == (1_ip + (kfl_paral_proc_node * ipart))) then
      !      IPARSLAVE = .TRUE.
      !      exit
      !   end if
      !end do

      if( IPARSLAVE ) then
!      if( INOTMASTER ) then
         icolor = 1_4
      else
         icolor = 2_4
      end if
      call PAR_COMM_SPLIT(icolor,PAR_COMM_SFC_WM,iproc,'IN MY CODE')
      PAR_MY_SFC_RANK_WM = int(iproc,ip)
      !
      ! Creates partitioners and master communicator
      !
      if( IPARSLAVE .or. IMASTER ) then
         icolor = 1_4
      else
         icolor = 2_4
      end if
 !     icolor = 1_4
      call PAR_COMM_SPLIT(icolor,PAR_COMM_SFC,iproc,'IN MY CODE')

      !
      ! Send sizes and stuff the slaves don't have
      !
      call PAR_BROADCAST(nelem,'IN SFC PARTITION WITH MASTER')
      call PAR_BROADCAST(npoin,'IN SFC PARTITION WITH MASTER')
      call PAR_BROADCAST(nelty,ngaus,'IN SFC PARTITION WITH MASTER')
#ifndef PNODE_VALUE
      call PAR_BROADCAST(mnode,'IN SFC PARTITION WITH MASTER')
#endif
      call PAR_COMM_RANK_AND_SIZE(PAR_COMM_SFC_WM,comm_rank,comm_size)
      nboxc=comm_size
      call PAR_MAX(nboxc,'IN SFC PARTITION WITH MASTER')

   end subroutine par_slaves_communicator

   !
   ! In this subroutine, the slaves will exchange the coordinates that they need in order to calculate their elements center of mass
   !
   subroutine par_exchange_coordinates

      use def_master,  only : IPARSLAVE
      use def_domain,  only : nnode

      implicit none

      integer(ip)             :: ielem,inode,ielty, ipoin, jpoin, indexi,aux_index
      integer(ip)             :: ipart, cpart
      integer(ip)             :: min_ipoin_needed, max_ipoin_needed
      integer(ip)             :: buffer_size_send, buffer_size_recv
      integer(ip)             :: npoin_local              ! Number of nodes after the exchange of needed coordinates

      integer(ip)             :: number_coord_need        ! total number of coordinates needed
      integer(ip), pointer    :: communication_table(:,:) ! table of communications, it says how many points are needed among processes
      integer(ip), pointer    :: max_coord_needed(:,:)    ! Max number of nodes that i need from the each slaves and the slave that has it
      integer(ip), pointer    :: ranking(:)               ! Ranking arrary to help create local numbering
      integer(ip), pointer    :: diction(:)               ! Dictionary to go from local to global numbering diction(new) = old
      type(i1p),   pointer    :: point_list_recv(:)       ! List with the points that i recv
      type(i1p),   pointer    :: point_list_send(:)       ! List with the points that i send
      type(r2p),   pointer    :: recv_coord_buffer(:)     ! Receive coordinates buffer
      type(r2p),   pointer    :: send_coord_buffer(:)     ! Send coordinates buffer

      real(rp),    pointer    :: dummr_buf2(:,:)          ! dummi real buffer (two entries)
      real(rp),    pointer    :: temp_coord(:,:)          ! temporal coordinates array

      logical(lg), pointer    :: bool_coord_needed(:)     ! booleans used to don't repeat searches

      !
      ! Pointers nullify
      !
      nullify(communication_table)
      nullify(max_coord_needed)
      nullify(ranking)

      nullify(point_list_recv)
      nullify(point_list_send)
      nullify(recv_coord_buffer)
      nullify(send_coord_buffer)

      nullify(dummr_buf2)
      nullify(temp_coord)
      nullify(bool_coord_needed)

      !
      ! Global pointers
      !
      nullify(lnods_loc)
      nullify(diction)
      nullify(invdict)

#ifndef MPI_OFF
      if( IPARSLAVE )then
         !
         ! Communication table allocation
         !
         call memory_alloca(par_memor,'communication_table','par_exchange_coordinates_sfc',communication_table,&
            & nboxc,nboxc,'DO_NOT_INITIALIZE',0_ip,0_ip)

         do ipart = 0_ip, nboxc-1_ip
            do cpart = 0_ip, nboxc-1_ip
               communication_table(ipart,cpart) = 0_ip
            end do
         end do
         !
         ! Check de list of nodes of my elements and determine if i have the coordinates of the node or i'll have to ask for them
         ! First pass around elements, check sizes
         !
         number_coord_need = 0_ip
         min_ipoin_needed  = npoin
         max_ipoin_needed  = 0_ip
         do ielem = 1, nelem_part

            ielty = abs(ltype_par(ielem))
            do inode =  1_ip, nnode(ielty)

               ipoin =  lnods_par(inode,ielem)
               if( ipoin < ipoi0 + 1_ip .or. ipoin > ipoi0 + npoin_part )then
                  !
                  ! The node is in not my range
                  !
                  number_coord_need = number_coord_need + 1_ip

                  if( ipoin < min_ipoin_needed ) min_ipoin_needed = ipoin
                  if( ipoin > max_ipoin_needed ) max_ipoin_needed = ipoin

               end if

            end do

         end do
         !
         ! Allocate boolean array to avoid repetition of coordinates
         !
         call memory_alloca(par_memor,' bool_coord_needed ',' par_exchange_coordinates_sfc ',&
            bool_coord_needed,max_ipoin_needed-min_ipoin_needed+1_ip)
         call memory_alloca(par_memor,' invdict           ',' par_exchange_coordinates_sfc ',&
            invdict,max( max_ipoin_needed, ipoi0+npoin_part)-min( min_ipoin_needed,ipoi0 ) + 1_ip)

         bool_coord_needed = .true.
         invdict           = -2500_ip
         !
         ! Allocate temporary array with the maximum number of coordinates that will be needed (normally oversized)
         !
         call memory_alloca(par_memor,' max_coord_needed ',' par_exchange_coordinates_sfc ',max_coord_needed,2_ip,mnode*nelem_part)
         !
         ! Second pass, mark the points i need and figure out who has them
         !
         number_coord_need = 0_ip

         do ielem = 1, nelem_part

            ielty = abs(ltype_par(ielem))
            do inode =  1_ip, nnode(ielty)

               ipoin =  lnods_par(inode,ielem)
               if( ipoin < ipoi0 + 1_ip .or. ipoin > ipoi0 + npoin_part )then
                  !
                  ! The node is in not my range
                  !
                  ! print*,  "DEBUG: soy ", PAR_MY_SFC_RANK_WM, "y yo no tengo el punto ", ipoin," que necesito :'(, sospecho que ", ipoin, " lo tiene ", ipoin/npoin_aver + min(1_ip, mod(ipoin,npoin_aver))-1_ip
                  if( bool_coord_needed(ipoin - min_ipoin_needed + 1_ip) )then

                     bool_coord_needed(ipoin - min_ipoin_needed + 1_ip) = .false.                                         ! won't count this point again

                     ipart = min(ipoin/npoin_aver + min(1_ip, mod(ipoin,npoin_aver))-1_ip,nboxc-1_ip)

                     number_coord_need = number_coord_need + 1_ip
                     max_coord_needed(1_ip,number_coord_need) = ipoin                                                     ! The ipoin i need
                     max_coord_needed(2_ip,number_coord_need) = ipart                                                     ! The rank of slave (partitioners communicator)


                     communication_table(PAR_MY_SFC_RANK_WM,ipart) = &
                        communication_table(PAR_MY_SFC_RANK_WM,ipart) + 1_ip


                  end if

               end if

            end do

         end do

         call memory_deallo(par_memor,' bool_coord_needed ',' par_exchange_coordinates_sfc ',bool_coord_needed)
         !
         ! Allocate a temporal array for the coordinates i will receive and store those i already have, also allocate a dictionary to create a local numbering
         !
         npoin_local = npoin_part+number_coord_need
         call memory_alloca(par_memor,' temp_coord ',' par_exchange_coordinates_sfc ',temp_coord,ndime,npoin_local)
         call memory_alloca(par_memor,' diction    ',' par_exchange_coordinates_sfc ',diction,npoin_local)
         do ipoin = 1_ip, npoin_part
            temp_coord(:,ipoin) = coord_par(:,ipoin)
            diction(ipoin)      = ipoi0 + ipoin
         end do
         aux_index = npoin_part  ! auxiliar index to store the coordinates i will receive
         !
         ! Reduction to assemble the communication table, partititoners + master
         !
         call PAR_SUM(communication_table, 'IN SFC PARTITION',INCLUDE_ROOT=.true.)

         call PAR_START_NON_BLOCKING_COMM(1_ip,nboxc-1_ip)
         !
         ! Communications loop over the partitioner slaves to say how many coordinates i need from them
         !
         call memory_alloca(par_memor,' point_list_send '  ,' par_exchange_coordinates_sfc ',point_list_send,  nboxc,'INITIALIZE',0_ip)
         call memory_alloca(par_memor,' point_list_recv '  ,' par_exchange_coordinates_sfc ',point_list_recv,  nboxc,'INITIALIZE',0_ip)
         call memory_alloca(par_memor,' recv_coord_buffer ',' par_exchange_coordinates_sfc ',recv_coord_buffer,nboxc,'INITIALIZE',0_ip)
         call memory_alloca(par_memor,' send_coord_buffer ',' par_exchange_coordinates_sfc ',send_coord_buffer,nboxc,'INITIALIZE',0_ip)
         do cpart = 0_ip, nboxc-1_ip

            if( PAR_MY_SFC_RANK_WM /= cpart )then
               !
               ! Prepare lists of points i need and i send
               !
               buffer_size_send = communication_table(PAR_MY_SFC_RANK_WM,cpart)
               buffer_size_recv = communication_table(cpart, PAR_MY_SFC_RANK_WM)

               call memory_alloca(par_memor,' point_list_send '  ,' par_exchange_coordinates_sfc ',point_list_send(cpart) % l, buffer_size_send)
               call memory_alloca(par_memor,' point_list_recv '  ,' par_exchange_coordinates_sfc ',point_list_recv(cpart) % l, buffer_size_recv)
               call memory_alloca(par_memor,' recv_coord_buffer ',' par_exchange_coordinates_sfc ',recv_coord_buffer(cpart) % a, ndime, buffer_size_send)
               call memory_alloca(par_memor,' send_coord_buffer ',' par_exchange_coordinates_sfc ',send_coord_buffer(cpart) % a, ndime, buffer_size_recv)

               indexi = 0_ip

               do ipoin = 1_ip, number_coord_need

                  if( max_coord_needed(2_ip,ipoin) == cpart )then
                     indexi = indexi + 1_ip
                     point_list_send(cpart) % l(indexi) = max_coord_needed(1_ip,ipoin)
                  end if

               end do

               if( indexi /= buffer_size_send ) then
                  print*,'Error with indexi=',indexi, buffer_size_send
                  call runend(' par_exchange_coordinates_sfc: size mismatch in points list exchange  ')
               end if
               !
               ! Order the points to send in order to be more efficient in future search
               !
               if( buffer_size_send /= 0_ip )call heapsorti1(2_ip, buffer_size_send, point_list_send(cpart) % l)
               !
               ! Send the list of points i need and receive the list my neighbour needs
               !
               call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
               call PAR_SEND_RECEIVE(point_list_send(cpart) % l ,point_list_recv(cpart) % l,'IN SFC PARTITION', cpart, 'NON BLOCKING')
               do ipoin = 1_ip, buffer_size_recv
                  ! print*,"DEBUG: ", PAR_MY_SFC_RANK_WM, cpart, point_list_recv(ipoin)
               end do

            end if

         end do
         call PAR_END_NON_BLOCKING_COMM(1_ip)

         call PAR_START_NON_BLOCKING_COMM(1_ip,nboxc-1_ip)
         do cpart = 0_ip, nboxc-1_ip

            if( PAR_MY_SFC_RANK_WM /= cpart )then
               !
               ! Prepare the list of coordinates the other slave needs, they are already ordered
               !
               buffer_size_recv = communication_table(cpart, PAR_MY_SFC_RANK_WM)
               indexi = 1_ip
               if( buffer_size_recv /= 0_ip )then

                  do ipoin = 1_ip, npoin_part
                     ! print*, "DEBUG: PAR_MY_SFC_RANK_WM ", PAR_MY_SFC_RANK_WM, cpart, " && ", ipoi0+ipoin," == ? " ,point_list_recv(indexi), " indexi ", indexi
                     if( ipoi0+ipoin == point_list_recv(cpart) % l(indexi) )then
                        send_coord_buffer(cpart) % a(:,indexi) = coord_par(:,ipoin)
                        indexi = indexi + 1_ip
                     end if
                     if( indexi == buffer_size_recv + 1_ip ) exit
                  end do

                  if( indexi-1_ip /= buffer_size_recv )then
                     print*, "RUNEND : PAR_MY_SFC_RANK_WM ", PAR_MY_SFC_RANK_WM, " indexi ", indexi-1_ip, " /= ", buffer_size_recv
                     call runend(' par_exchange_coordinates_sfc: size mismatch in coordinates exchange  ')
                  end if

               end if
               !
               ! Send and receive coordinates
               !
               call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
               call PAR_SEND_RECEIVE(send_coord_buffer(cpart) % a, recv_coord_buffer(cpart) % a, 'IN SFC PARTITION', cpart, 'NON BLOCKING')

               do ipoin = 1_ip,buffer_size_send
                  ! write(99000+PAR_MY_SFC_RANK_WM,*) "DEBUG: PAR_MY_SFC_RANK_WM ", PAR_MY_SFC_RANK_WM, " ipoin ", ipoin, point_list_send(cpart) % l(ipoin), " coordinates i have ", recv_coord_buffer(cpart) % a (:,ipoin)
               end do

            end if

         end do  ! Communications loop
         call PAR_END_NON_BLOCKING_COMM(1_ip)
         !
         ! Store the received coords
         !
         indexi     = 1_ip
         aux_index = npoin_part
         do cpart = 0_ip, nboxc-1_ip

            do ipoin = 1_ip, communication_table(PAR_MY_SFC_RANK_WM,cpart) ! buffer_size_send
               temp_coord(:,aux_index + ipoin) = recv_coord_buffer(cpart) % a(:,ipoin)
               diction(aux_index + ipoin)      = point_list_send(cpart) % l(ipoin)
            end do
            aux_index  = aux_index + communication_table(PAR_MY_SFC_RANK_WM,cpart)

            if( PAR_MY_SFC_RANK_WM /= cpart )then
               indexi = indexi + communication_table(PAR_MY_SFC_RANK_WM,cpart)
            else
               indexi = indexi + npoin_part

            end if

         end do
         !
         ! Set a local numbering to search more efficiently among coordinates and create the local coordinates array
         !
         call memory_alloca(par_memor,' ranking '  ,' par_exchange_coordinates_sfc ',ranking,npoin_local)
         call memory_alloca(par_memor,' coord_loc ',' par_exchange_coordinates_sfc ',coord_loc,ndime,npoin_local)
         call memory_alloca(par_memor,' lnods_loc ',' par_exchange_coordinates_sfc ',lnods_loc,mnode,nelem_part)

         !
         ! initialization of the ranking array, this will tell you the ranks in the dictionary
         !
         do ipoin = 1_ip, npoin_local
            ranking(ipoin) = ipoin
         end do

         do ipoin = 1_ip, npoin_local
            ! write(99000+PAR_MY_SFC_RANK_WM,*) "DEBUG A: PAR_MY_SFC_RANK_WM ", PAR_MY_SFC_RANK_WM, " total coordinates ", ipoin, diction(ipoin), ranking(ipoin), temp_coord(:,ipoin)
         end do
         !
         ! Order the dictionary and the ranking array
         !
         call heapsorti2(2_ip,npoin_local,diction,ranking)

         !
         ! Write your local coordinates (finally!)
         !
         do ipoin = 1_ip, npoin_local

            coord_loc(:,ipoin) = temp_coord(:,ranking(ipoin))
            jpoin              = diction(ipoin) - min(min_ipoin_needed,ipoi0) + 1_ip
            !print*, "DEBUG:  PAR_MY_SFC_RANK_WM invdict ", PAR_MY_SFC_RANK_WM, ipoin,jpoin,min_ipoin_needed
            invdict(jpoin)     = ipoin
         end do
         do ipoin = 1_ip, npoin_local
            ! print*, "DEBUG:  PAR_MY_SFC_RANK_WM invdict ", PAR_MY_SFC_RANK_WM, ipoin, invdict(ipoin)
         end do
         !
         ! Write lnods local
         !
         !print*, "DEBUG: searching ...", PAR_MY_SFC_RANK_WM

         do ielem = 1_ip, nelem_part

            ielty = abs(ltype_par(ielem))
            do inode =  1_ip, nnode(ielty)

               ipoin =  lnods_par(inode,ielem)
               lnods_loc(inode, ielem) = invdict(ipoin - min(min_ipoin_needed,ipoi0) + 1_ip)
               ! print*, "DEBUG:  PAR_MY_SFC_RANK_WM invdict ", PAR_MY_SFC_RANK_WM, ipoin, invdict(ipoin)
               !              !
               !              ! Determine the owner of this point
               !              !
               !              ipart = ipoin/npoin_aver + min(1_ip, mod(ipoin,npoin_aver))-1_ip
               ! print*, "DEBUG: ", index, aux_index, ipart
               ! call array_search(diction, npoin_local, ipoin, jpoin)
               ! call array_search(diction(index:aux_index), aux_index-index+1_ip, ipoin, jpoin)
               ! print*, "DEBUG: PAR_MY_SFC_RANK_WM ",PAR_MY_SFC_RANK_WM, " global ipoin ", ipoin,jpoin

               ! lnods_loc(inode, ielem) = jpoin
               ! bool_coord_needed(ipoin) = .false.

            end do

         end do

         call memory_deallo(par_memor,' diction '          ,' par_exchange_coordinates_sfc ',diction          )
         call memory_deallo(par_memor,' point_list_send '  ,' par_exchange_coordinates_sfc ',point_list_send  )
         call memory_deallo(par_memor,' point_list_recv '  ,' par_exchange_coordinates_sfc ',point_list_recv  )
         call memory_deallo(par_memor,' recv_coord_buffer ',' par_exchange_coordinates_sfc ',recv_coord_buffer)
         call memory_deallo(par_memor,' send_coord_buffer ',' par_exchange_coordinates_sfc ',send_coord_buffer)
         call memory_deallo(par_memor,' bool_coord_needed ',' par_exchange_coordinates_sfc ',bool_coord_needed)

         if(associated(lnods_par)) call memory_deallo(par_memor,' lnods_par ',' par_exchange_coordinates_sfc ',lnods_par)
         if(associated(communication_table)) call memory_deallo(par_memor,' communication_table ',&
            'par_exchange_coordinates_sfc ',communication_table)
         if(associated(max_coord_needed)) call memory_deallo(par_memor,' max_coord_needed ',' max_coord_needed ',max_coord_needed)
         if(associated(ranking)) call memory_deallo(par_memor,' ranking ',' par_exchange_coordinates_sfc ',ranking)
         if(associated(temp_coord)) call memory_deallo(par_memor,' temp_coord ',&
            ' par_exchange_coordinates_sfc',temp_coord)

         call PAR_BARRIER('IN SFC PARTITION')

      end if
#endif

   end subroutine par_exchange_coordinates
   !
   ! In this subroutine, master will send coords, nodes, elements and basic information to partition the mesh in parallel
   ! to a subgroup of partitioners slaves simulating parallel reading of the mesh
   !
   subroutine par_distribute_mesh()

      use def_master,  only : IPARSLAVE

      implicit none

      integer(ip)             :: ielem,ipoin,ipart,indexi
      integer(ip)             :: message_size, inisend, finsend
      integer(ip)             :: nelem_aver      ! Number of elements of all partitions except last one
      integer(ip)             :: dummi(1_ip)
      integer(ip), pointer    :: dummi_buf2(:,:) ! dummi integer buffer (two entries)
      real(rp),    pointer    :: dummr_buf2(:,:) ! dummi real buffer (two entries)
      integer(ip), pointer    :: send_lnods_buff(:,:)
      real(rp),    pointer    :: send_coord_buff(:,:) ! sending coord buffer

#ifndef MPI_OFF
      nullify(coord_par)
      nullify(ltype_par)
      nullify(lnods_par)
      nullify(send_coord_buff)
      nullify(send_lnods_buff)

      nullify(dummi_buf2)
      nullify(dummr_buf2)

      nullify(nelem_part_master4)
      nullify(npoin_part_master4)
      !
      ! Determine the number of elements and points that each partitioner slave will treat
      !
      !if( mod(npoin, nboxc ) >= nboxc / 2_ip )then
      !   nelem_aver = nelem / nboxc + 1_ip             ! Number of elements of all partitions except last one
      !   npoin_aver = npoin / nboxc + 1_ip             ! Number of points of all partitions except last one
      !else
         nelem_aver = nelem / nboxc
         npoin_aver = npoin / nboxc
      !end if


      if( PAR_MY_SFC_RANK_WM == nboxc - 1_ip ) then
         !
         ! I am the last salve and have a different load
         !
         nelem_part = nelem - nelem_aver*(nboxc - 1_ip)
         npoin_part = npoin - npoin_aver*(nboxc - 1_ip)
         !print*, "DEBUG:nelem_aver PAR_MY_SFC_RANK_WM ",PAR_MY_SFC_RANK_WM, nelem_part, npoin_part, nelem, npoin
      else
         nelem_part = nelem_aver
         npoin_part = npoin_aver
         !print*, "DEBUG:nelem_aver PAR_MY_SFC_RANK_WM ",PAR_MY_SFC_RANK_WM, nelem_part, npoin_part, nelem, npoin
      end if

      !
      ! Master sends the messages for the slaves.
      !
      if( IMASTER )then

         call memory_alloca(par_memor,' nelem_part_master ',' par_distribute_mesh ' ,nelem_part_master4,int(nboxc+1_ip,4_ip),'DO_NOT_INITIALIZE',0_4)
         call memory_alloca(par_memor,' npoin_part_master ',' par_distribute_mesh ' ,npoin_part_master4,int(nboxc+1_ip,4_ip),'DO_NOT_INITIALIZE',0_4)
         nelem_part_master4 = 0_4
         npoin_part_master4 = 0_4

         do ipart = 1_ip, nboxc-1_ip

            inisend      = ( ipart-1_ip ) * nelem_aver + 1_ip
            finsend      = ipart * nelem_aver
            !
            ! ltype
            !
            message_size = nelem_aver
            nelem_part_master4(ipart) = int(message_size,4_ip)
            call PAR_SEND_RECEIVE(message_size,0_ip, ltype(inisend:finsend), dummi, 'IN SFC PARTITION WITH MASTER', ipart)
            !
            ! lnods
            !
            call memory_alloca(par_memor,' send_lnods_buff ',' par_distribute_mesh ' ,send_lnods_buff,mnode,message_size)
            indexi = 1_ip
            do ielem = inisend, finsend
               send_lnods_buff(:,indexi) = lnods(:,ielem)
               indexi = indexi + 1_ip
            end do
            call PAR_SEND_RECEIVE(send_lnods_buff, dummi_buf2, 'IN SFC PARTITION WITH MASTER', ipart)
            call memory_deallo(par_memor,' send_lnods_buff ',' par_distribute_mesh ' ,send_lnods_buff)
            !
            ! coordinates
            !
            inisend      = ( ipart-1_ip ) * npoin_aver + 1_ip
            finsend      = ipart * npoin_aver
            message_size = finsend-inisend+1_ip
            npoin_part_master4(ipart) = int(message_size,4_ip)
            call memory_alloca(par_memor,' send_coord_buff ',' par_distribute_mesh ' ,send_coord_buff,ndime,message_size)
            indexi = 1_ip
            do ipoin = inisend, finsend
               send_coord_buff(:,indexi) = coord(:,ipoin)
               indexi = indexi + 1_ip
            end do
            call PAR_SEND_RECEIVE(send_coord_buff, dummr_buf2, 'IN SFC PARTITION WITH MASTER', ipart)
            call memory_deallo(par_memor,' send_coord_buff ',' par_distribute_mesh ' ,send_coord_buff)

         end do
         !
         ! For the last slave
         !
         inisend      = ( nboxc-1_ip ) * nelem_aver + 1_ip
         finsend      = nelem
         !
         ! ltype
         !
         message_size = finsend - inisend + 1_ip ! nelem_aver +  ( nelem - nelem_aver*(nboxc - 1_ip) )
         nelem_part_master4(nboxc) = int(message_size,4_ip)
         call PAR_SEND_RECEIVE(message_size,0_ip,ltype(inisend:finsend), dummi, 'IN SFC PARTITION WITH MASTER', nboxc)
         !
         ! lnods
         !
         call memory_alloca(par_memor,' send_lnods_buff ',' par_distribute_mesh ' ,send_lnods_buff,mnode,message_size)
         indexi = 1_ip
         do ielem = inisend, finsend
            send_lnods_buff(:,indexi) = lnods(:,ielem)
            indexi = indexi + 1_ip
         end do
         call PAR_SEND_RECEIVE(send_lnods_buff, dummi_buf2, 'IN SFC PARTITION WITH MASTER', ipart)
         !
         ! coordinates
         !
         inisend      = ( nboxc-1_ip ) * npoin_aver + 1_ip
         finsend      = npoin

         message_size = finsend - inisend+1_ip
         npoin_part_master4(nboxc) = int(message_size,4_ip)
         call memory_alloca(par_memor,' send_coord_buff ',' par_distribute_mesh ' ,send_coord_buff,ndime,message_size)
         indexi = 1_ip
         do ipoin = inisend, finsend
            send_coord_buff(:,indexi) = coord(:,ipoin)
            indexi = indexi + 1_ip
         end do
         call PAR_SEND_RECEIVE(send_coord_buff, dummr_buf2, 'IN SFC PARTITION WITH MASTER', ipart)

      else if( IPARSLAVE )then

         !
         ! Setting starting indexes to go from my local ipoin/ielem to the global one, initially PAR_MY_CODE_RANK is used for simplicity,
         ! but the rank in the partition communicator has to be used instead
         !
         ipoi0 = ( PAR_MY_CODE_RANK - 1_ip ) * npoin_aver
         !
         ! Partitioner slaves allocate memory to receive the messages
         !
         call memory_alloca(par_memor,'coord_par','par_distribute_mesh',coord_par,ndime,npoin_part)
         call memory_alloca(par_memor,'ltype_par','par_distribute_mesh',ltype_par,nelem_part)
         call memory_alloca(par_memor,'lnods_par','par_distribute_mesh',lnods_par,mnode,nelem_part)
         !
         ! The slaves receive the messages
         ! ltype
         !
         call PAR_SEND_RECEIVE(0_ip, nelem_part, dummi, ltype_par, 'IN SFC PARTITION WITH MASTER', 0_ip)
         !
         ! lnods
         !
         call PAR_SEND_RECEIVE(dummi_buf2, lnods_par, 'IN SFC PARTITION WITH MASTER', 0_ip)
         !
         ! coordinates
         !
         call PAR_SEND_RECEIVE(dummr_buf2, coord_par, 'IN SFC PARTITION WITH MASTER', 0_ip)
      end if
      ! print*, "DEBUG: PAR_MY_CODE_RANK ", PAR_MY_CODE_RANK, " ltype_par ", ltype_par
      if( IMASTER )then
         do ipoin = 1_ip, npoin
            ! write(98000+PAR_MY_CODE_RANK,*) "DEBUG: PAR_MY_CODE_RANK ", PAR_MY_CODE_RANK, " coord ", coord(:,ipoin), " ipoin ", ipoin, npoin
         end do
      end if
      if( IPARSLAVE )then
         do ielem = 1_ip, nelem_part
            ! write(97000+PAR_MY_CODE_RANK,*) "DEBUG: PAR_MY_CODE_RANK ", " lnods_par ", lnods_par(:,ielem), " ielem ", ielem
         end do
         do ipoin = 1_ip, npoin_part
            ! write(97000+PAR_MY_CODE_RANK,*) "DEBUG: PAR_MY_CODE_RANK ", PAR_MY_CODE_RANK," coord_par ", coord_par(:,ipoin), " ipoin ", ipoin, ipoi0,npoin_part
         end do
      end if
#endif


   end subroutine par_distribute_mesh
   !
   ! Generate the SFC partition input
   !
   subroutine par_sfc_input()

      character(100), PARAMETER :: vacal = "par_sfc_input"
      real(rp)                  :: cent_mass(3)
      integer(ip)               :: ienti
      integer(ip)               :: ielty, inode, ipoin, iweig

      if( IPARSLAVE ) then

         !
         ! Allocate memory
         !
         if(sfc_criteria >= 1_ip) then
            nenti=nelem_part
         else
            nenti=npoin_part
         endif

         call memory_alloca(par_memor,' lenti ',vacal,lenti,ndime,nenti)
         call memory_alloca(par_memor,' lweig ',vacal,lweig,nenti)

         if(sfc_criteria >= 1_ip) then

            do ienti = 1_ip, nenti
               !
               ! Evaluate mass center (mc)
               !
               cent_mass(1:3) = 0.0_rp
               ielty = abs(ltype_par(ienti))
               do inode = 1_ip, nnode(ielty)
                  ipoin = lnods_loc(inode,ienti)
                  cent_mass(1:ndime) = cent_mass(1:ndime) + coord_loc(1:ndime,ipoin)
               end do
               cent_mass(1:ndime) = cent_mass(1:ndime) / real(nnode(ielty),rp)
               lenti(1:ndime,ienti) = cent_mass(1:ndime)
               !
               ! evaluate weight
               !
               if(sfc_criteria == 2_ip) then
                  iweig = ngaus(abs(ltype_par(ienti)))
               else
                  iweig = 1_rp
               endif
               lweig(ienti) = iweig
            end do
         else

            do ienti = 1_ip, nenti
               lenti(1:ndime,ienti) = coord_par(1:ndime,ienti)
               lweig(ienti) = 1_ip
            end do

         endif


      endif

   end subroutine par_sfc_input

   subroutine par_deallocate

      character(100), PARAMETER :: vacal = "par_dealocate_sfc"

      if(associated(coord_loc)) call memory_deallo(par_memor,' coord_loc ',vacal,coord_loc)
      if(associated(coord_par)) call memory_deallo(par_memor,' coord_par ',vacal,coord_par)
      if(associated(ltype_par)) call memory_deallo(par_memor,' ltype_par ',vacal,ltype_par)
      if(associated(lnods_loc)) call memory_deallo(par_memor,' lnods_loc ',vacal,lnods_loc)

   end subroutine par_deallocate
   !
   !
   !
   subroutine gather_to_master(lepar_par,npart)

      use mod_communications,  only : PAR_GATHERV

      implicit none

      integer(ip), pointer, intent(inout)  :: lepar_par(:)
      integer(ip),          intent(in)     :: npart
      integer(ip), pointer                 :: dummy(:)
      integer(ip)                          :: ielem
      integer(ip)                          :: nmaxp

      character(100), PARAMETER :: vacal = "gather_to_master"

      integer(ip), allocatable          :: resul(:)
      nmaxp = npart


      nullify(dummy)
      call memory_alloca(par_memor,' dummy ',vacal,dummy,1_ip)

      if( IPARSLAVE )then

            call PAR_GATHERV(lepar_par,dummy,nelem_part_master4,'IN SFC PARTITION WITH MASTER')

      end if

      if( IMASTER )then
         nullify(lepar_par)
         call memory_alloca(par_memor,' lepar_par ',vacal,lepar_par,nelem)
         lepar_par  = -1_ip
         if(sfc_criteria >= 1) then
            call PAR_GATHERV(dummy,lepar_par,nelem_part_master4,'IN SFC PARTITION WITH MASTER')

         else if(sfc_criteria == 0 ) then
            call runend("PAR_SFC: option not uploaded yet")
            !Master must assign a partition to the elements according to the nodes partition
         endif
         !
         !  CHECK:
         !

         do ielem = 1, size(lepar_par)
            if(lepar_par(ielem)<0 .or. lepar_par(ielem) > npart ) then
               call runend("Bad partition")
            endif
         enddo
         if(sfc_check == 1) then
            allocate(resul(npart))
            resul(:)=0_ip
            do ielem = 1, size(lepar_par)
               resul(lepar_par(ielem))=resul(lepar_par(ielem))+1_ip
            enddo
            print *, "Result: ", resul
            print *, "Lepar:  ", lepar_par
            print *, "        "
            deallocate(resul)
            stop
         endif
      end if
      call memory_deallo(par_memor,' dummy ',vacal,dummy)

   end subroutine gather_to_master
   !
   ! Partition visualization
   !
   subroutine par_visualization_sfc(lepar,visop,dboxc,dim_bin_core)

      use def_domain,              only : coord,ltype,lnods,nnode
      use mod_space_filling_curve, only : space_filling_curve_output
      use mod_maths,               only : maths_mapping_coord_to_3d

      implicit none

      integer(ip), pointer, intent(in) :: lepar(:)
      integer(ip),          intent(in) :: visop     ! visualization option 0: no visualize, 1: weights, 2: partition
      integer(ip),          intent(in) :: dboxc(3)  ! #bin boxes per direction in the coarse bin (coherent with partition)
      integer(ip),          intent(in) :: dim_bin_core ! dimendion of local grid used for partition (coherent with parition)
      integer(ip), pointer             :: visfd(:)
      integer(ip)                      :: ielem,iboxe,idime,inode,ielty,ipoin
      integer(ip)                      :: bcoof(3), nboxf
      real(rp)                         :: min_coord(3)    ! bounding box min
      real(rp)                         :: max_coord(3)    ! bounding bos max
      real(rp)                         :: cent_mass(3)
      integer(ip)                      :: dboxf(3)
      real(rp)                         :: deltf(3)        ! edge sizes of the fine bin boxes

      character(100), PARAMETER :: vacal = "par_visualization_sfc"

      if(visop > 0_ip .and. visop > 2_ip) then
         call runend("par_visualization_sfc:invalid visop option")
      endif

      if( IMASTER .and. visop > 0_ip ) then
         do idime = 1_ip,ndime
            max_coord(idime) = maxval(coord(idime,:))
            min_coord(idime) = minval(coord(idime,:))
         end do

         dboxf(1:3) = dim_bin_core*dboxc(1:3)
         nboxf = dboxf(1)**ndime

         deltf(3)=1.0_rp
         do idime = 1_ip, ndime
            deltf(idime) = (max_coord(idime) - min_coord(idime)) / real(dboxf(idime),rp)
         end do

         nullify(visfd)
         call memory_alloca(par_memor,' visfd ',vacal,visfd,nboxf)

         do ielem = 1_ip, nelem
            cent_mass(1:3) = 0.0_rp
            ielty = abs(ltype(ielem))
            do inode = 1_ip,nnode(ielty)
               ipoin = lnods(inode,ielem)
               cent_mass(1:ndime) = cent_mass(1:ndime) + coord(1:ndime,ipoin)
            end do
            cent_mass(1:ndime) = cent_mass(1:ndime) / real(nnode(ielty),rp)

            call maths_mapping_coord_to_3d(ndime,dboxf,min_coord,max_coord,cent_mass,bcoof(1),bcoof(2),bcoof(3))
            bcoof(1:ndime) = min(bcoof(1:ndime),dboxf(1:ndime))
            iboxe = dboxf(2_ip)*dboxf(1_ip)*(bcoof(3)-1_ip) + dboxf(1_ip)*(bcoof(2)-1_ip) + bcoof(1)

            if(visop == 1_ip) then
               visfd(iboxe) = lepar(ielem)
            else if (visop == 2_ip) then
               visfd(iboxe) = visfd(iboxe) + 1_ip
            endif
         end do
         call space_filling_curve_output(ndime,dboxf(1),deltf(1_ip),deltf(2_ip),deltf(3_ip),visfd,'LEXICAL')
         call memory_deallo(par_memor,' visfd ',vacal,visfd)

      end if


   end subroutine par_visualization_sfc


end module mod_redistribute
