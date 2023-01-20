!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Repartitioning
!> @{
!> @file    mod_repartitioning.f90
!> @author  houzeaux
!> @date    2019-06-18
!> @brief   Module for repartitioning
!> @details Mesh repartitioning
!-----------------------------------------------------------------------

module mod_repartitioning

   use def_master
   use def_kermod
   use def_parall
   use mod_messages
   use mod_memory
   use mod_communications
   use mod_parall,            only : PAR_COMM_MY_CODE
   use mod_par_partitioning,  only : par_partitioning
   use mod_messages,          only : messages_live
   use def_parall,            only : kfl_weigh_par
   use mod_parall,            only : par_memor
   use mod_parall,            only : lun_repar_par
   use mod_parall,            only : PAR_WEIGHT_OFF       
   use mod_parall_destructor, only : parall_destructor  
   use def_master,            only : INOTMASTER
   use mod_mpio_par_postpr,   only : init_redistrib_post
   use mod_mpio_par_postpr,   only : init_redistrib_rst
   use mod_iofile,            only : iofile_flush_unit
   use mod_communications,    only : PAR_ALLGATHER
   use mod_parall,            only : PAR_CODE_SIZE
   use mod_parall,            only : PAR_MY_CODE_RANK
   use mod_outfor,            only : outfor
   use mod_output,            only : output_open_files
   use mod_maths_basic,       only : maths_weighted_linear_regression
   implicit none
   private

   character(100), PARAMETER :: vacal = "mod_repartitioning"
   !
   ! Repartitionig
   !
   integer(ip)       :: num_passes = 0
   real(rp)          :: time_ave   = 0.0_rp
   real(rp)          :: time_max   = 0.0_rp
   real(rp)          :: my_time    = 0.0_rp
   real(rp)          :: init_time  = 0.0_rp
   real(rp)          :: time_destructor
   real(rp)          :: time_partitioning
   real(rp)          :: time_redistribution
   real(rp)          :: time_domain
   real(rp)          :: time_total
   real(rp), pointer :: lcorr(:)

   ! 
   ! SFC partition corrections
   !
   integer(ip)              :: nunit = 0_ip
   real(rp),    pointer     :: lcorr_ant(:)
   real(rp),    pointer     :: reg_loads(:)
   real(rp),    pointer     :: reg_times(:)
   real(rp),    pointer     :: reg_weigh(:)
   
   ! Input paramters
   integer(ip)              :: method
   integer(ip)              :: window
   real(rp)                 :: minload 
   real(rp)                 :: maxload 
   real(rp)                 :: toler
   real(rp)                 :: wfact

   public :: repartitioning
   public :: repartitioning_initialization
   public :: repartitioning_timings

contains

   !-----------------------------------------------------------------------
   !> 
   !> @author  houzeaux
   !> @date    2019-10-28
   !> @brief   Output timings
   !> @details Output timings of repartitioning
   !> 
   !-----------------------------------------------------------------------

   subroutine repartitioning_timings()

      ioutp(1) = num_repart_par
      routp(1) = time_total
      routp(2) = time_destructor
      routp(3) = time_partitioning
      routp(4) = time_domain
      routp(5) = time_redistribution

      call outfor(103_ip)

   end subroutine repartitioning_timings

   !-----------------------------------------------------------------------
   !> 
   !> @author  houzeaux
   !> @date    2019-06-18
   !> @brief   Initialization
   !> @details Initialization of repartitioning variables
   !> 
   !-----------------------------------------------------------------------

   subroutine repartitioning_initialization

      num_repart_par      = 0
      num_passes          = 0
      nunit               = 0
      time_ave            = 0.0_rp
      time_max            = 0.0_rp
      my_time             = 0.0_rp
      time_destructor     = 0.0_rp
      time_partitioning   = 0.0_rp
      time_redistribution = 0.0_rp
      time_domain         = 0.0_rp
      time_total          = 0.0_rp
      nullify(lcorr)
      nullify(lcorr_ant)

   end subroutine repartitioning_initialization

   !-----------------------------------------------------------------------
   !> 
   !> @author  houzeaux
   !> @date    2019-06-18
   !> @brief   Main subroutine for repartitioning
   !> @details Repartition the mesh and redistribute data
   !> 
   !-----------------------------------------------------------------------
   subroutine repartitioning() 

      integer(ip)         :: imodu
      real(rp)            :: lb
      real(rp)            :: time1,time2,timei,timef

      if( kfl_gotim == 0 ) return
      num_passes = num_passes + 1
      nullify(lcorr)
      do_repartitioning = .false.

      if( kfl_repart_par /= 0 ) then

         call cputim(timei)

         !-----------------------------------------------------------------
         !
         ! Compute load balance
         !
         !-----------------------------------------------------------------

         if( kfl_repart_module_par > 0 .and. kfl_repart_criterion_par /= -1 ) then 
            imodu        = kfl_repart_module_par
            my_time      = my_time + cpu_modul(kfl_repart_criterion_par,imodu) - init_time
            time_ave     = my_time
            time_max     = my_time
            call PAR_MAX(time_max)
            call PAR_AVERAGE(time_ave)
            lb = time_ave / time_max

            if( IMASTER ) then
               if( lb < repart_threshold_lb .and. num_passes == kfl_repart_freq_par ) then
                  write(lun_repar_par,'(i7,3(1x,e13.6),1x,i1)') num_passes,time_ave,time_max,lb,1_ip
               else
                  write(lun_repar_par,'(i7,3(1x,e13.6),1x,i1)') num_passes,time_ave,time_max,lb,0_ip             
               end if
               call iofile_flush_unit(lun_repar_par)
            end if
         else
            imodu = 0
            lb    = -1.0_rp
         end if

         !-----------------------------------------------------------------
         !
         ! Redistribute
         !
         !-----------------------------------------------------------------

         if( lb < repart_threshold_lb .and. num_passes == kfl_repart_freq_par ) then
            num_repart_par    = num_repart_par + 1
            !kfl_weigh_par     = PAR_WEIGHT_OFF
            do_repartitioning = .true.
            call messages_live('---')
            call messages_live('REPARTITIONING','START SECTION')
            if( lb > 0.0_rp ) call messages_live('LOAD BALANCE= '//trim(retost(lb,REAL_FORMAT='(f6.2)')))
            call memory_allocation_mode('DEALLOCATE BEFORE ALLOCATING')
            !
            ! Destructors
            !
            call cputim(time1)
            call messages_live('DESTRUCT SOLVERS')
            call solver_destructor()
            call messages_live('DESTRUCT DOMAIN')
            call coupli_destructor()
            call domain_destructor()
            call parall_destructor()
            call output_open_files(REPARTITIONING=.true.)
            call cputim(time2)
            time_destructor = time_destructor + time2 - time1
            time1 = time2
            !
            ! Compute weights and repartition
            !
            if( kfl_repart_criterion_par /= -1 ) then
               call repartitioning_weights(my_time, lcorr,lb)          
               call par_partitioning(READ_MESH=.false.,LCORRECT=lcorr)
               call repartitioning_weights_deallocate(lcorr)
            else        
               call par_partitioning(READ_MESH=.false.)
            end if
            call cputim(time2)
            time_partitioning = time_partitioning + time2 - time1
            time1 = time2
            !
            ! Construct domain 
            !
            init_redistrib_post = .false.
            init_redistrib_rst  = .false.
            call Domain()
            call cputim(time2)
            time_domain = time_domain + time2 - time1
            time1 = time2
            call Redist()             ! Variables to be redistributed
            call cputim(time2)
            time_redistribution = time_redistribution + time2 - time1
            time1 = time2
            call Solmem()             ! Variables to be reallocated
            call Begrun()
            !
            ! Reinitialize
            !
            call memory_allocation_mode('DO NOT DEALLOCATE BEFORE ALLOCATING')
            call messages_live('REPARTITIONING','END SECTION')
            call messages_live('---')
            if( kfl_repart_module_par > 0 .and. kfl_repart_criterion_par /= -1 ) then
               init_time  = cpu_modul(kfl_repart_criterion_par,imodu) 
               my_time    = 0.0_rp
            end if
            num_passes = 0

         end if

         call cputim(timef)
         time_total = time_total + timef - timei

      end if

   end subroutine repartitioning

   !-----------------------------------------------------------------------
   !> 
   !> @author  borrell
   !> @date    2019-06-18
   !> @brief   Evaluation of weights for partition correction
   !> @details Evaluation of weights for partition correction
   !> 
   !-----------------------------------------------------------------------
   subroutine repartitioning_weights(eltim, lcorr,lb)

      use mod_debug  !rick: to disapear

      real(rp),             intent(in)    :: eltim       ! elaptsed time
      real(rp),   pointer,  intent(inout) :: lcorr(:)    ! corrections per subdomain
      real(rp),             intent(in)    :: lb          ! rick: to disapear
      real(rp)                            :: sumti,allti,aveti
      real(rp)                            :: a,b
      real(rp)                            :: sumlcorr, sumlcorr2
      real(rp)                            :: alpha
      real(rp),   pointer                 :: auxr(:), lelap(:)
      integer(ip)                         :: ii,ind0,ind1
      integer(ip)                         :: iprev, inext
      integer(ip)                         :: P,rank
      integer(ip)                         :: usereg
      integer(ip), pointer                :: fixed(:)   

      usereg=0.0
      if( INOTMASTER ) then
         
         if( kfl_repart_par == 1 ) then

            !
            ! Initialization
            !
            P      = PAR_CODE_SIZE-1
            rank   = PAR_MY_CODE_RANK
            if( .not. associated(lcorr)) call memory_alloca(par_memor,'lcorr',vacal,lcorr,P)
            if(nunit == 0) then
               call repartitioning_weights_initialization(rank)
               call memory_alloca(par_memor,'lcorr_ant',vacal,lcorr_ant,P)
               lcorr_ant = 1.0_rp
            endif

            !
            ! Get all elapsed times
            !
            nullify(lelap,fixed)
            call memory_alloca(par_memor,'fixed',vacal,fixed,P)
            call memory_alloca(par_memor,'lelap',vacal,lelap,P)
            call PAR_ALLGATHER(eltim,lelap,'IN MY CODE WITHOUT MASTER')
            allti = sum(lelap)
            if(allti<=0.0_rp) call runend('repartitioning_weights: the accumulated elapsed time can not be zero')
            aveti = allti/real(P,rp)
   
            !
            ! Evaluate new coefficients
            !
            if( method == CPU_REPART_LOAD ) then
               

               do ii=1,P
                  if(abs(lelap(ii)/aveti-1.0)>toler) then
                     lcorr(ii) = lcorr_ant(ii)*aveti/lelap(ii)
                     lcorr(ii) = max(lcorr(ii),minload)
                     lcorr(ii) = min(lcorr(ii),maxload)
                     fixed(ii) = 0_ip
                  else
                     fixed(ii) = 1_ip
                     lcorr(ii) = lcorr_ant(ii)
                  endif
               enddo

            elseif( method == CPU_REPART_OLR .or. method == CPU_REPART_WLR) then
               !
               ! Offsets within the window
               !
               ind0 = mod(nunit-1,window)+1
               ind1 = mod(nunit,window)+1
               fixed(:) = 0_ip
               ! 
               ! Evaluate new relative accumulated time and regresson weigh
               !
               sumti = 0.0_rp
               do ii = 1, rank
                  sumti = sumti + lelap(ii)
               enddo
               reg_times(ind0)=sumti/aveti
               if(method == CPU_REPART_OLR) then
                  reg_weigh(ind0)=1.0 
               else
                  sumti = 0.0_rp
                  reg_weigh(ind1)=wfact
                  iprev = ind1
                  do ii=1, window-1
                     inext = mod(iprev,window)+1
                     reg_weigh(inext)=wfact*reg_weigh(iprev) 
                     iprev = inext   
                  enddo
               endif
               !
               ! Linear regression
               !
               call maths_weighted_linear_regression(reg_loads,reg_times,nunit,a,b,reg_weigh,window)
               !
               ! Intersection
               !
               if( abs(a) > epsilon(1.0_rp) ) then
                  reg_loads(ind1) =  (real(rank,rp)-b)/a
               else
                  reg_loads(ind1) =  reg_loads(ind0)
               endif
               !
               ! Check crosses
               !
               nullify(auxr)
               call memory_alloca(par_memor,'auxr',vacal,auxr,P)
               call PAR_ALLGATHER(reg_loads(ind1),auxr,'IN MY CODE WITHOUT MASTER')
               usereg = 1_ip
               loop_ii : do ii=2,P
                  if( auxr(ii) <= auxr(ii-1) ) then
                     usereg = 0_ip
                     exit loop_ii
                  endif
               enddo loop_ii
               !
               ! Evaluate new lcorr coeficients
               !
               !if(nunit < window/2+1) usereg = 0_ip  !use first load update rick check
               if(usereg > 0 ) then 
                  !
                  ! Use linear regretion
                  !
                  lcorr(1)=auxr(1)
                  do ii=2,P
                     lcorr(ii) = auxr(ii) - auxr(ii-1)
                     lcorr(ii) = max(lcorr(ii),minload)
                     lcorr(ii) = min(lcorr(ii),maxload)
                  enddo
               else
                  !
                  ! Reescale weight proportionaly
                  !
                  do ii=1,P
                     lcorr(ii) = lcorr_ant(ii)*aveti/lelap(ii)
                     lcorr(ii) = max(lcorr(ii),minload)
                     lcorr(ii) = min(lcorr(ii),maxload)
                  enddo
               endif
               call memory_deallo(par_memor,'auxr',vacal,auxr)
            endif
            !
            ! Normalization respeting minload and maxload
            !
            sumlcorr = 0.0_rp
            sumlcorr2 = 0.0_rp
            do ii=1,P
               if(fixed(ii)==0) then
                  sumlcorr = sumlcorr + lcorr(ii)
               else
                  sumlcorr2 = sumlcorr2 + lcorr(ii)
               endif
            enddo

            if( sum(fixed) /= P ) then
               alpha    = (P-sumlcorr2) / sumlcorr
               do ii=1,P
                  if(fixed(ii)==0) then
                     lcorr(ii) = alpha * lcorr(ii)
                     lcorr(ii) = max(lcorr(ii),minload)
                     lcorr(ii) = min(lcorr(ii),maxload)
                  endif
               enddo
               !
               ! Final normalization
               !
               sumlcorr = 0.0_rp
               do ii=1,P
                  if(fixed(ii)==0) then
                     sumlcorr = sumlcorr + lcorr(ii)
                  endif
               enddo
               alpha    = (P-sumlcorr2) / sumlcorr
               do ii=1,P
                  if(fixed(ii)==0) then
                     lcorr(ii)=  alpha * lcorr(ii)
                     lcorr_ant(ii) = lcorr(ii)
                  endif
               enddo
            end if

            if( method == CPU_REPART_OLR .or. method == CPU_REPART_WLR) then
               !
               ! Evaluate reg_loads 
               !
               sumlcorr = 0
               do ii = 1, rank
                  sumlcorr = sumlcorr + lcorr(ii)
               enddo
               reg_loads(ind1)=sumlcorr
               nunit = nunit + 1_ip
            endif

            call memory_deallo(par_memor,'lelap',vacal,lelap)
            call memory_deallo(par_memor,'fixed',vacal,fixed)


         else if( kfl_repart_par == 2 ) then
            call runend("repartitioning_weights: kfl_repart_par==2 not encoded yet")

         end if
      end if
      call PAR_MAX_ALL(usereg)
      if(usereg>0) then
         call messages_live('WEIGHTS FROM LINEAR REGRESSION')
      else
         call messages_live('WEIGHTS RESCALED')
      endif

   end subroutine repartitioning_weights

   subroutine repartitioning_weights_initialization(rank)

      integer(ip), intent(in) :: rank

      nunit        = 2_ip
      window       = kfl_repart_windo_par   
      minload      = repart_min_weight
      maxload      = repart_max_weight
      toler        = repart_toler
      wfact        = repart_wfact
      method       = kfl_repart_method_par

      nullify(reg_loads,reg_times)
      call memory_alloca(par_memor,'reg_loads',vacal,reg_loads,window)
      call memory_alloca(par_memor,'reg_times',vacal,reg_times,window)
      call memory_alloca(par_memor,'reg_weigh',vacal,reg_weigh,window)

      reg_loads(1) = 0.0_rp
      reg_times(1) = 0.0_rp        
      reg_weigh(1) = 1.0_rp    !rick: veure que fem amb weights
      reg_loads(2) = real(rank,rp)              

   end subroutine repartitioning_weights_initialization 

   !-----------------------------------------------------------------------
   !> 
   !> @author  houzeaux
   !> @date    2019-09-02
   !> @brief   Deallocate
   !> @details Deallocate arrays of repartition
   !> 
   !-----------------------------------------------------------------------

   subroutine repartitioning_weights_deallocate(lcorr)

      real(rp),   pointer,  intent(inout)  :: lcorr(:)

      call memory_deallo(par_memor,'lcorr',vacal,lcorr)

   end subroutine repartitioning_weights_deallocate

end module mod_repartitioning
