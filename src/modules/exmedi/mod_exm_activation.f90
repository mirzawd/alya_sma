!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_exm_activation
   use def_kintyp_basic, only : i1p, ip, rp, lg 
   implicit none

   private

   integer(ip), parameter :: array_block_size = 50_ip ! the pts_stim_reach arrays are allocated in chunks of this number of elements

   !If a table of stimuli is used, pts_stim_reach(i) % l(j) variable will store the id of j-th point for i-th stimulus
   !Each slave has it's own arrays
   type(i1p),   dimension(:), pointer :: pts_stim_reach  => null() ! access with % l(i)
   integer(ip), dimension(:), pointer :: npts_stim_reach => null() ! number of points in each of the above subarrays

   public :: exm_stim_enumerate_activated_nodes ! Identify points ot be activated. To be called before timesteps start only once.
   public :: exm_stim_compute                   ! Calculate stimuli on every timestep
   public :: exm_stim_exchange                  ! Initial exchange of related variables, to put into sendat. Call once!
   public :: exm_stim_synchronize               ! Synchronize variables modified during the timesteps between slaves. Call on every timestep
   public :: exm_stim_allocate                  ! Allocate memory for stimulus table

   private :: exm_stim_add_point
   private :: exm_stim_allocate_reach_array
   private :: exm_stim_resize_container
   private :: exm_stim_compact_container
contains


subroutine exm_stim_allocate()
   use def_exmedi, only : aptim, apval_exm, aplap_exm, aprea_exm, apcen_exm, nstim_exm
   use mod_memory, only : memory_alloca
   use def_master, only : mem_modul, modul
   implicit none
   integer(ip)  :: nstim_exm_tmp

   nstim_exm_tmp = nstim_exm
   if (nstim_exm<0_ip) then 
      nstim_exm_tmp = 1_ip ! for BOUNDARY, TODO: TEST
   end if


   if (nstim_exm_tmp>0_ip) then
      call memory_alloca(mem_modul(1:2,modul),'apval_exm', 'mod_exm_activation', apval_exm, nstim_exm_tmp)
      call memory_alloca(mem_modul(1:2,modul),'aplap_exm', 'mod_exm_activation', aplap_exm, nstim_exm_tmp)
      call memory_alloca(mem_modul(1:2,modul),'aprea_exm', 'mod_exm_activation', aprea_exm, nstim_exm_tmp)
      call memory_alloca(mem_modul(1:2,modul),'aptim',     'mod_exm_activation', aptim,     nstim_exm_tmp)
      call memory_alloca(mem_modul(1:2,modul),'apcen_exm', 'mod_exm_activation', apcen_exm, 3_ip, nstim_exm_tmp)
      apval_exm = 0.0_rp
      aplap_exm = 0.0_rp
      aprea_exm = 0.0_rp
      aptim     = 0.0_rp
      apcen_exm = 0.0_rp
   end if

end subroutine exm_stim_allocate



subroutine exm_stim_exchange()
   use def_exmedi, only : aptim, apval_exm, aplap_exm, aprea_exm, apcen_exm, nstim_exm
   use mod_exchange, only : exchange_add, exchange_end, exchange_init
   implicit none

   if ( nstim_exm .ne. 0_ip ) then
      call exchange_init()
      call exchange_add(aptim)
      call exchange_add(apval_exm)
      call exchange_add(aplap_exm)
      call exchange_add(aprea_exm)
      call exchange_add(apcen_exm)
      call exchange_end()
   end if
end subroutine exm_stim_exchange

subroutine exm_stim_synchronize()
   use def_exmedi,         only : aptim, nstim_exm, kfl_appty_exm
   use mod_communications, only : PAR_MIN
   implicit none
   integer(ip) :: istim

   if (kfl_appty_exm == 2_ip) then ! if FLASH
      !aptim >= 0 normal stimulus set in exm.dat, not changed during the execution, 
      !With FLASH, aptim  is modified in the code, set to -1 to kill the stimulus. Needs synchronization
      !Therefore taking min to pull the changes from the slaves
      do istim=1,nstim_exm
         call PAR_MIN(aptim(istim))
      end do
   end if
end subroutine exm_stim_synchronize


subroutine exm_stim_enumerate_activated_nodes()
   use def_exmedi,   only : nstim_exm, apcen_exm, aprea_exm, modst_exm
   use def_domain,   only : npoin, coord, ndime
   use mod_messages, only : messages_live
   use def_master,   only : intost, INOTMASTER, retost
   implicit none

   integer(ip) :: istim, ipoin
   real(rp) :: sqreach
   real(rp) :: xicen, yicen, zicen, dicen
      
   if (modst_exm>=0_ip) then !only if field is not used for the activation
      call messages_live("EXMEDI LOCATES POINTS CLOSE TO THE ACTIVATION STIMULI CENTERS")
      
      if (INOTMASTER ) then 
         if ( nstim_exm > 0_ip ) then !BOUNDARY not included
            call exm_stim_allocate_reach_array(nstim_exm)

            do istim= 1, nstim_exm

               sqreach= aprea_exm(istim)*aprea_exm(istim)    
               do ipoin= 1,npoin
                  xicen = coord(1,ipoin)-apcen_exm(1,istim)
                  yicen = coord(2,ipoin)-apcen_exm(2,istim)
                  zicen = 0.0_rp
                  if (ndime == 3) zicen = coord(ndime,ipoin)-apcen_exm(ndime,istim)

                  if ( (abs(xicen)<=aprea_exm(istim)) .and. &
                     (abs(yicen)<=aprea_exm(istim)) .and.   &
                     (abs(zicen)<=aprea_exm(istim)) ) then

                     dicen = xicen*xicen + yicen*yicen + zicen*zicen

                     if (dicen <= sqreach) then
                        call exm_stim_add_point(istim, ipoin)
                     end if
                  end if
               end do

               call exm_stim_compact_container(istim)

               if ( npts_stim_reach(istim)==0_ip ) then
                  call messages_live("NO MESH NODES WITHIN "//trim(retost(aprea_exm(istim)))//" OF STIMULUS POINT "//trim(intost(istim)),'WARNING')
               end if
            end do
         end if
      end if
   end if
end subroutine exm_stim_enumerate_activated_nodes


subroutine exm_stim_add_point(istim, ipoin)
   use mod_memory, only : memory_size
   implicit none
   integer(ip), intent(in) :: istim, ipoin

   npts_stim_reach(istim) = npts_stim_reach(istim) + 1

   !if we reached the end of the array - reallocate
   call exm_stim_resize_container(istim, npts_stim_reach(istim))

   pts_stim_reach(istim) % l( npts_stim_reach(istim) ) = ipoin
end subroutine exm_stim_add_point


subroutine exm_stim_allocate_reach_array(nstim)
   use mod_memory, only : memory_alloca
   use def_master, only : mem_modul, modul
   implicit none
   
   integer(ip), intent(in) :: nstim

   call memory_alloca( mem_modul(1:2,modul),'pts_stim_reach', 'mod_exm_activation', pts_stim_reach,  nstim)
   call memory_alloca( mem_modul(1:2,modul),'npts_stim_reach', 'mod_exm_activation', npts_stim_reach,  nstim)
   
   npts_stim_reach = 0_ip 

end subroutine exm_stim_allocate_reach_array

subroutine exm_stim_resize_container(istim, nelements)
   use mod_memory, only : memory_alloca, memory_size, memory_deallo
   use def_master, only : mem_modul, modul, intost
   implicit none

   integer(ip), intent(in)             :: istim, nelements
   integer(ip)                         :: reach_size, new_size, i
   integer(ip), dimension(:), pointer  :: temp_array => null()
   logical(lg)                         :: array_exists = .false.

   if( istim > memory_size(npts_stim_reach) ) then
      call runend("exm_stim_resize_container: istim = "//trim(intost(istim))//" is outside allocated memory for "//trim(intost(memory_size(npts_stim_reach)))//" stimuli")
   end if

   reach_size = 0_ip

   if ( associated(pts_stim_reach(istim) % l) ) then 
      reach_size = memory_size(pts_stim_reach(istim) % l)
      array_exists = .true.
   end if

   if ( nelements>reach_size ) then
      new_size = ceiling(real(nelements)/real(array_block_size), kind=ip) * array_block_size !don't care about real precision

      if ( array_exists ) then 

         call memory_alloca( mem_modul(1:2,modul),'temp_array', 'mod_exm_activation', temp_array, reach_size )

         do i=1,reach_size
            temp_array(i) = pts_stim_reach(istim) % l(i)
         end do
         
         call memory_deallo( mem_modul(1:2,modul),'pts_stim_reach % l', 'mod_exm_activation', pts_stim_reach(istim) % l )
         call memory_alloca( mem_modul(1:2,modul),'pts_stim_reach % l', 'mod_exm_activation', pts_stim_reach(istim) % l,  new_size )
         pts_stim_reach(istim) % l(:) = 0_ip

         do i=1,reach_size
            pts_stim_reach(istim) % l(i) = temp_array(i)
         end do

         call memory_deallo( mem_modul(1:2,modul),'temp_array', 'mod_exm_activation', temp_array )

      else 

         call memory_alloca( mem_modul(1:2,modul),'pts_stim_reach % l', 'mod_exm_activation', pts_stim_reach(istim) % l,  new_size )
      end if
   end if 

end subroutine exm_stim_resize_container


subroutine exm_stim_compact_container(istim)
   use mod_memory, only : memory_alloca, memory_size, memory_deallo
   use def_master, only : mem_modul, modul, intost
   implicit none

   integer(ip), intent(in)             :: istim
   integer(ip)                         :: reach_size, i, nelements
   integer(ip), dimension(:), pointer  :: temp_array => null()

   if( istim > memory_size(npts_stim_reach) ) then
      call runend("exm_stim_compact_container: istim = "//trim(intost(istim))//" is outside allocated memory for "//trim(intost(memory_size(npts_stim_reach)))//" stimuli")
   end if


   if ( associated(pts_stim_reach(istim) % l) ) then 
      reach_size = memory_size(pts_stim_reach(istim) % l)
      nelements  = npts_stim_reach(istim)

      if ( nelements .ne. reach_size ) then
         call memory_alloca( mem_modul(1:2,modul),'temp_array', 'mod_exm_activation', temp_array, nelements )

         do i=1,nelements
            temp_array(i) = pts_stim_reach(istim) % l(i)
         end do
         
         call memory_deallo( mem_modul(1:2,modul),'pts_stim_reach % l', 'mod_exm_activation', pts_stim_reach(istim) % l )
         call memory_alloca( mem_modul(1:2,modul),'pts_stim_reach % l', 'mod_exm_activation', pts_stim_reach(istim) % l,  nelements )

         do i=1,nelements
            pts_stim_reach(istim) % l(i) = temp_array(i)
         end do

         call memory_deallo( mem_modul(1:2,modul),'temp_array', 'mod_exm_activation', temp_array )
      end if 
   end if

end subroutine exm_stim_compact_container





!------------------------------------------------------------------------
! Compute of stimulus current
! Setup the appfi_exm (applied current field) to the EP problem
! The parameters are being read in exm_reaphy.f90 
!-----------------------------------------------------------------------
! TODO: refactor this function
subroutine exm_stim_compute()

   use      def_master
   use      def_domain
   use      def_exmedi
   implicit none

   integer(ip) :: ipoin,istim,iboun,pblty,pnodb,inodb
   logical(lg) :: flash = .false.
   logical(lg) :: istim_compute = .true.
 
   integer(ip) :: nstim_exm_tmp

   call times(1) % ini()
   
   if (INOTMASTER) then
      flash = ( kfl_appty_exm == 2_ip )


      !  ltime = cutim -dtime                        !!! <<<---  Iapp is that of the CURRENT time step
      appfi_exm = 0.0_rp

      if (modst_exm < 0_ip) then
         !
         ! Stimuli in a field
         !
         do ipoin= 1,npoin
            if (xfiel(-modst_exm) % a(2,ipoin,1) .ge. 0.0_rp) then
               if (cutim .ge. xfiel(-modst_exm) % a(2,ipoin,1)) then              
                  appfi_exm(ipoin) = xfiel(-modst_exm) % a(1,ipoin,1)
                  if (cutim .gt. (xfiel (-modst_exm) % a (2,ipoin,1) + xfiel (-modst_exm) % a (3,ipoin,1))) then 
                     appfi_exm(ipoin) = 0.0_rp
                  end if
               end if
            end if
         end do

      else  ! not field and there is an activation table
         nstim_exm_tmp = nstim_exm
         if (nstim_exm<0_ip) then ! BOUNDARY, TODO:TEST
            nstim_exm_tmp = 1_ip
         end if

         do istim= 1, nstim_exm_tmp

            istim_compute = .false.
            if (aptim(istim) >= 0_ip) then
               if (cutim >= aptim(istim) .and. cutim <= (aplap_exm(istim)+ aptim(istim))) then
                  istim_compute = .true.

                  if (flash) then !if flash, kill the stimulus for the next time step
                     aptim(istim) = -HUGE(1.0_rp)
                  end if
               end if
            end if

            if (istim_compute) then 
                  
               qneto_exm = 0.0_rp ! This resets QNET on every beat      

               if (nstim_exm < 0_ip) then !BOUNDARY
                  do iboun= 1,nboun
                     pblty = ltypb(iboun) 
                     pnodb = nnode(pblty)
      
                     if (lbset(iboun)== nstis_exm) then
      
                        do inodb= 1,pnodb
                           ipoin= lnodb(inodb,iboun)
      
                           if (cutim >= aptim(istim) .and. cutim <= (aplap_exm(istim)+ aptim(istim))) then
                              appfi_exm(ipoin) = apval_exm(1)
                              if (flash) then
                                 elmag(ipoin,1)= apval_exm(1)      
                                 elmag(ipoin,2)= apval_exm(1)      
                                 elmag(ipoin,3)= apval_exm(1)
                              end if
                           end if
      
                        end do
                     end if
                  end do
               else !not BOUNDARY
                  do ipoin= 1,npts_stim_reach(istim)
                     appfi_exm( pts_stim_reach(istim) % l(ipoin) ) = apval_exm(istim)  
                     if (flash) then !flash is always on volatge
                        elmag    (pts_stim_reach(istim) % l(ipoin),1_ip) = apval_exm(istim)
                        elmag    (pts_stim_reach(istim) % l(ipoin),2_ip) = apval_exm(istim)
                        elmag    (pts_stim_reach(istim) % l(ipoin),3_ip) = apval_exm(istim)
                     end if
                  end do   
               end if ! (nstim_exm < 0_ip)

            end if !istim_compute
         end do !istim

      end if !activation table

   end if

   call times(1) % add()

end subroutine exm_stim_compute


end module mod_exm_activation
