!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Output
!> @{
!> @file    mod_outfor.f90
!> @author  houzeaux
!> @date    2018-12-14
!> @brief   Output
!> @details Module for Alya output. All messages are grouped here
!>          in order to easily change Alya format
!-----------------------------------------------------------------------

module mod_outfor

#include "def_vector_size.inc"
  use def_kintyp,            only : ip,rp
  use mod_communications,    only : PAR_AVERAGE
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_SUM
  use mod_communications,    only : PAR_MIN
  use mod_memory,            only : memory_unit
  use mod_maths,             only : maths_unit
  use mod_element_data_base, only : element_data_base_memory
  use mod_messages,          only : livinf
  use mod_messages,          only : messages_live
  use mod_elmgeo,            only : element_type
  use mod_iofile,            only : iofile_flush_unit
  use def_master,            only : coutp
  use def_master,            only : ioutp
  use def_master,            only : routp
  use def_master,            only : namod
  use def_master,            only : READ_AND_RUN
  use def_master,            only : INOTSLAVE
  use def_master,            only : INOTMASTER
  use def_master,            only : ISLAVE
  use def_master,            only : ISEQUEN
  use def_master,            only : IMASTER
  use def_master,            only : iblok
  use def_master,            only : dtime
  use def_master,            only : ittim
  use def_master,            only : itcou
  use def_master,            only : modul
  use def_master,            only : itinn
  use def_master,            only : kfl_timco
  use def_master,            only : kfl_async
  use def_master,            only : kfl_rstar
  use def_master,            only : zeror
  use def_master,            only : npart
  use def_master,            only : title
  use mod_strings,           only : integer_to_string
  use mod_direct_solver
  use mod_parall
  use def_parall
  use def_solver
  use def_domain
  use def_kermod
  use mod_log
  use mod_run_config,        only : run_config
  use mod_live_info_config,  only : live_info_config
  
  implicit none

  public :: outfor

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-14
  !> @brief   Output messages in files
  !> @details Output messages in Alya files
  !> 
  !-----------------------------------------------------------------------

  subroutine outfor(itask,nunit,messa,REAL_NUMBERS,INT_NUMBERS,INT_LIST,CHA_LIST)

    integer(ip),                          intent(in) :: itask                                   !< What to output
    integer(ip),                optional, intent(in) :: nunit                                   !< Unit where to output
    character(len=*),           optional, intent(in) :: messa                                   !< An optional message
    real(rp),                   optional, intent(in) :: REAL_NUMBERS(:)                         !< Optional list of real numbers
    integer(ip),                optional, intent(in) :: INT_NUMBERS(:)                          !< Optional list of integers
    integer(ip),  pointer,      optional, intent(in) :: INT_LIST(:)                             !< Optional pointer for list of integers
    character(len=*),           optional, intent(in) :: CHA_LIST(:)                             !< Optional pointer for list of integers
    integer(ip)                                      :: jtask,ii,nsolv,isolv,memor_value4,dummi
    integer(ip)                                      :: dumm1,dumm2
    integer(ip)                                      :: ngpu
    integer(ip),  save                               :: ipass=0
    real(rp),     save                               :: xtime=0.0_rp
    real(rp)                                         :: dummr
    real(rp)                                         :: unit_factor,unit_dot_factor
    character(20)                                    :: me010,elem1,elem2
    character(100)                                   :: me100
    character(40)                                    :: woutl='--------------------------------------'
    character(3)                                     :: wmpi3,wblas
    character(1)                                     :: unit_char,unit_dot_char
    integer(8)                                       :: memor_value
    character(6)                                     :: memor_char
    real(rp)                                         :: memor_factor
    real(rp)                                         :: mem_nzmat_min,mem_nzmat_max,mem_nzmat_ave
    real(rp)                                         :: flops_max,flops_ave
    real(rp)                                         :: flops_dot_max,flops_dot_ave
    integer(ip)                                      :: value_mix(4)
    integer(ip)                                      :: value_ave(2)
    integer(ip)                                      :: nzmat_min,nzmat_max,nzmat_ave  
    integer(ip)                                      :: omp_chunk_size_min             
    integer(ip)                                      :: omp_chunk_size_max             
    integer(ip)                                      :: omp_chunk_size_ave 
    integer(ip)                                      :: omp_scheduling(4)
    integer(4)                                       :: istat4
    integer(4)                                       :: nunit4
    integer(4)                                       :: lun_outpu4

    if( INOTSLAVE .or. itask <= 0 ) then
       !
       ! What to do
       !
       jtask      = abs(itask)
       if( present(nunit) ) then
          nunit4  = int(nunit,4)
       else
          nunit4  = int(lun_outpu,4)
       end if
       lun_outpu4 = int(lun_outpu,4)
       !
       ! Optional arguments
       !
       if( present(REAL_NUMBERS) ) then
          ii = size(REAL_NUMBERS)
          routp(1:ii) = REAL_NUMBERS(1:ii)
       end if

       if( present(INT_NUMBERS) ) then
          ii = size(INT_NUMBERS)
          ioutp(1:ii) = INT_NUMBERS(1:ii)
       end if

       select case (jtask)

       case(0_ip)
         call log_runend_error(messa)

       case(1_ip)
          !
          ! Error message
          !
          if( present(messa) ) then
             if( nunit4 == 0 ) then
                write(lun_outpu4,100) trim(messa)
                call iofile_flush_unit(lun_outpu4)
             else
                write(nunit4,100) trim(messa)
                call iofile_flush_unit(nunit4)
             end if
          end if

       case(2_ip)
          !
          ! Warning message
          !
          if( present(messa) ) then
             if( nunit4 == 0 ) then
                write(lun_outpu4,200) trim(messa)
             else
                call messages_live(trim(messa),'WARNING')
                write(nunit4,200)     trim(messa)
             end if
          end if

       case(3_ip)
          !
          ! Warning in module MODUL
          !
          write(lun_outpu4,300)   'WARNINGS HAVE BEEN FOUND IN '//&
               trim(namod(modul))//' DATA FILE'
          write(lun_outpu4,*)
          call iofile_flush_unit(lun_outpu4)
          write(nunit4,*)
          call messages_live('WARNINGS HAVE BEEN FOUND IN MODULE '//trim(namod(modul)))

       case(4_ip)
          !
          ! Error in module MODUL
          !
          if( present(messa) ) then
             if( trim(messa) == '1' ) then
                call runend(trim(messa)//' ERROR HAS BEEN FOUND IN '//&
                     trim(namod(modul))//' DATA FILE')
             else
                call runend(trim(messa)//' ERRORS HAVE BEEN FOUND IN '//&
                     trim(namod(modul))//' DATA FILE')
             end if
          end if

       case(5_ip)
          !
          ! Begin of an inner iteration (Solite)
          !
          if( solve_sol(1) % kfl_solve == 1 ) then
             write(solve_sol(1) % lun_solve,500) ittim,itcou,itinn(modul)
             call iofile_flush_unit(solve_sol(1) % lun_solve)
             return
          end if

       case(6_ip)
          !
          ! End of run (Turnof)
          !
          write(nunit4,600) trim(namod(modul))

       case(7_ip)
          !
          ! Time step (Begste)
          !
          !if(run_config%log) write(nunit4,700) ittim

       case(8_ip)
          !
          ! Critical time step (*_tistep)
          !
          if(run_config%log) then
             if( kfl_timco/=2 ) then
                if( ioutp(1)/=0.and.ioutp(2)/=1 ) then
                   if( routp(1)>0.0_rp ) then
                      write(nunit4,800) namod(modul),routp(1),dtime/routp(1)
                   else
                      write(nunit4,800) namod(modul),routp(1),0.0_rp
                   end if
                end if
             else
                if( ioutp(1)/=0.and.ioutp(2)/=1)&
                     write(nunit4,801) namod(modul),routp(1),routp(2),routp(3)
             end if
          end if

       case(9_ip)
          !
          ! Coupling convergence (*_concou)
          !
          if(run_config%log) write(nunit4,900) coutp(1),routp(1)

       case(10_ip)
          !
          ! Coupling convergence (Concou)
          !
          if(run_config%log) write(nunit4,1000) itcou,iblok,routp(1),routp(1)-xtime
          xtime = routp(1)

       case(11_ip)
          !
          ! Restart run (*_outinf)
          !
          write(nunit4,1100)

       case(12_ip)
          !
          ! Title (*_outinf)
          !
          if( kfl_rstar == 2 ) then
             write(UNIT=nunit4,FMT=1100)
          else
             write(UNIT=nunit4,FMT=1200) namod(modul),trim(title)
          end if

       case(13_ip)
          !
          ! Wrong type of element
          !
          me010=integer_to_string(int(nunit4,ip))
          call runend('A WRONG ELEMENT WITH '//trim(me010)//' NODES HAS BEEN FOUND')

       case(14_ip)
          !
          ! End of analysis (runend)
          !
          call log_end_of_analysis()

       case(15_ip)
          !
          ! Time step (setgts)
          !
          if(run_config%log) write(nunit4,1500) ioutp(1),routp(1),routp(2)

       case(16_ip)
          !
          ! header (openfi)
          !
          call runend("USE MOD_WRITE_FILES SUBROUTINES")

       case(17_ip)
          !
          ! Information (outinf)
          !
          !write(nunit4,'(a)') '* Problem data'
          if( coutp(1) == '') coutp(1)='NONE'
          elem1 = integer_to_string(ioutp(3))
          elem2 = integer_to_string(ioutp(4))
          write(nunit4,1700)&
               trim(coutp(1)),trim(coutp(2)),trim(coutp(3)),&
               routp(5),routp(6),routp(7),&
               ioutp(1),&
               bandw_dom,profi_dom,&
               routp(1),routp(2),routp(3),trim(elem1),&
               routp(4),trim(elem2),ioutp(5),ioutp(6),ioutp(7),&
               trim(coutp(4)),trim(coutp(5)),ioutp(8),&
               ioutp(9),ioutp(10),&
               trim(coutp(6)),trim(coutp(7)),&
               ioutp(2),ioutp(11)
          if( kfl_edge_elements == 1 ) then
             write(nunit4,1701) meshe(ndivi) % nedge
          end if
          write(nunit4,*)

       case(18_ip)
          !
          ! CPU times (outcpu)
          !
          !write(nunit4,'(a)') '* Computing times'
          write(nunit4,1800)&
               routp( 1),&                ! Total time
               routp( 2),routp( 3),&      ! Starting operations
               routp( 4),routp( 5),&      ! Geometry reading
               routp( 6),routp( 7),&      ! Set reading
               routp( 8),routp( 9),&      ! Boundary condition reading
               routp(20),routp(21),&      ! Fields reading
               routp(10),routp(11),&      ! Mesh partitionining
               routp(14),routp(15),&      ! Mesh multiplication
               routp(16),routp(17),&      ! Domain construction
               &   routp(22),routp(23),&
               &   routp(24),routp(25),&
               &   routp(26),routp(27),&
               &   routp(28),routp(29),&
               &   routp(30),routp(31),&
               routp(18),routp(19)        ! Additional arrays

       case(19_ip)
          !
          ! CPU times (outcpu)
          !
          write(nunit4,1900)&
               trim(coutp(1)),&
               routp( 1) , routp( 2)     ! Total time

       case(20_ip)
          !
          ! CPU times (outcpu)
          !
          !write(nunit4,2000)&
          !     routp(1),routp(2),routp(3),routp(4)
          write(nunit4,2001)&
               routp(1),routp(2)

       case(21_ip)
          !
          ! Memory (outmem)
          !
          !if( ipass1 == 0 ) write(nunit4,'(a)') '* Memory'
          !ipass1 = ipass1 + 1
          if( ioutp(50) == 1 ) then
             write(nunit4,2100)
             write(nunit4,2102)
          else
             write(nunit4,2103)
          end if
          write(nunit4,2101) &
               routp(1),trim(coutp(1)),&
               routp(2),trim(coutp(2)),&
               routp(3),trim(coutp(3))

       case(22_ip)
          !
          ! Memory (outmem)
          !
          write(nunit4,2200) trim(coutp(1)),routp(1),trim(coutp(2))

       case(23_ip)
          !
          ! Memory (outmem)
          !
          return

       case(24_ip)
          !
          ! Memory (outmem)
          !
          write(nunit4,2400) routp(1),trim(coutp(1))

       case(25_ip)
          !
          ! Generic section titles for modules log files
          !
          if( present(messa) ) then
             me100 = ''
             dummi = len(trim(messa),kind=ip)
             do ii=1,min(dummi,len(me100,kind=ip))
                me100(ii:ii)='-'
             end do
             write(nunit4,2500) trim(me100),trim(messa),trim(me100)
          end if

       case(26_ip)
          !
          ! End of run (Turnof)
          !
          return

       case(27_ip)
          !
          ! Title (*_outinf)
          !
          return

       case(28_ip)
          !
          ! Steady state (*_cvgunk)
          !
          write(nunit4,2800) namod(modul),ittim

       case(29_ip)
          !
          ! Summary of computing time: title
          !
          write(nunit4,2900) max(routp(1),0.0_rp)

       case(30_ip)
          !
          ! Summary of computing time: section
          !
          write(nunit4,3000) adjustl(coutp(1)),max(routp(1),0.0_rp),max(routp(2),0.0_rp)

       case(31_ip)
          !
          ! Summary of computing time: subsection
          !
          write(nunit4,3100) adjustl(coutp(1)),max(routp(1),0.0_rp),max(routp(2),0.0_rp)

       case(32_ip)
          !
          ! Sets results: iteration header
          !
          write(nunit4,3200) ioutp(1),ioutp(2),ioutp(3),routp(1)

       case(33_ip)
          !
          ! Sets results: number of columns
          !
          write(nunit4,3300) ioutp(1),ioutp(2)

       case(34_ip)
          !
          ! Set results: file header
          !
          write(nunit4,3400) trim(coutp(1))

       case(35_ip)
          !
          ! Set results: variable, column
          !
          write(nunit4,3500) trim(coutp(1)),ioutp(1)

       case(36_ip)
          !
          ! Solver header
          !
          write(nunit4,3600)

       case(37_ip)
          !
          ! Solver statistics
          !
          nsolv=size(solve_sol,kind=ip)
          if( nsolv >= 1 .and. INOTSLAVE ) write(nunit4,3600)
          do isolv = 1,nsolv
             if( solve_sol(isolv) % nsolv > 0 ) then
                if( solve_sol(isolv) % kfl_algso/=-999 ) then
                   !
                   ! SpMV statistics
                   ! CPU_SPMV(1) = Computations
                   ! CPU_SPMV(2) = Communications
                   ! CPU_SPMV(3) = Total
                   !
                   solve_sol(isolv) % cpu_spmv(3) = solve_sol(isolv) % cpu_spmv(1) + solve_sol(isolv) % cpu_spmv(2)

                   solve_sol(isolv) % cpu_spmv(4) = solve_sol(isolv) % cpu_spmv(1) ! Comput.
                   solve_sol(isolv) % cpu_spmv(5) = solve_sol(isolv) % cpu_spmv(1) ! Comput.
                   solve_sol(isolv) % cpu_spmv(6) = solve_sol(isolv) % cpu_spmv(2) ! Comm.
                   solve_sol(isolv) % cpu_spmv(7) = solve_sol(isolv) % cpu_spmv(2) ! Comm.
                   solve_sol(isolv) % cpu_spmv(8) = solve_sol(isolv) % cpu_spmv(3) ! Total
                   solve_sol(isolv) % cpu_spmv(9) = solve_sol(isolv) % cpu_spmv(3) ! Total
                   !
                   ! Flops of SpMV
                   ! One multiple and add (*,+) per coefficient
                   !
                   flops_ave = 0.0_rp
                   flops_max = 0.0_rp
                   if( INOTMASTER .and. solve_sol(isolv) % num_spmv > 0 ) then
                      flops_max = 2.0_rp * real(solve_sol(isolv) % nzmat,rp) &
                           * real(solve_sol(isolv) % num_spmv,rp) &
                           / ( solve_sol(isolv) % cpu_spmv(3) + zeror )
                      flops_ave = flops_max
                   end if

                   call PAR_AVERAGE(solve_sol(isolv) % cpu_spmv(4))
                   call PAR_AVERAGE(solve_sol(isolv) % cpu_spmv(6))
                   call PAR_AVERAGE(solve_sol(isolv) % cpu_spmv(8))
                   call PAR_AVERAGE(flops_ave)
                   call PAR_MAX    (solve_sol(isolv) % cpu_spmv(5))
                   call PAR_MAX    (solve_sol(isolv) % cpu_spmv(7))
                   call PAR_MAX    (solve_sol(isolv) % cpu_spmv(9))
                   call PAR_MAX    (flops_max)
                   !
                   ! Dot product statistics
                   !
                   solve_sol(isolv) % cpu_dot(3) = solve_sol(isolv) % cpu_dot(1) + solve_sol(isolv) % cpu_dot(2)

                   solve_sol(isolv) % cpu_dot(4) = solve_sol(isolv) % cpu_dot(1) ! Comput.
                   solve_sol(isolv) % cpu_dot(5) = solve_sol(isolv) % cpu_dot(1) ! Comput.
                   solve_sol(isolv) % cpu_dot(6) = solve_sol(isolv) % cpu_dot(2) ! Comm.
                   solve_sol(isolv) % cpu_dot(7) = solve_sol(isolv) % cpu_dot(2) ! Comm.
                   solve_sol(isolv) % cpu_dot(8) = solve_sol(isolv) % cpu_dot(3) ! Total
                   solve_sol(isolv) % cpu_dot(9) = solve_sol(isolv) % cpu_dot(3) ! Total
                   !
                   ! Flops of dot product
                   ! One multiple and add (*,+) per unknown
                   !
                   flops_dot_ave = 0.0_rp
                   flops_dot_max = 0.0_rp
                   if( INOTMASTER .and. solve_sol(isolv) % num_dot > 0 ) then
                      flops_dot_max = 2.0_rp * real(solve_sol(isolv) % nequa_own*solve_sol(isolv) % ndofn,rp) &
                           * real(solve_sol(isolv) % num_dot,rp) &
                           / ( solve_sol(isolv) % cpu_dot(3) + zeror )
                      flops_dot_ave = flops_dot_max
                   end if

                   !call runend('O.K.!')
                   call PAR_AVERAGE(solve_sol(isolv) % cpu_dot(4))
                   call PAR_AVERAGE(solve_sol(isolv) % cpu_dot(6))
                   call PAR_AVERAGE(solve_sol(isolv) % cpu_dot(8))
                   call PAR_AVERAGE(flops_dot_ave)
                   call PAR_MAX    (solve_sol(isolv) % cpu_dot(5))
                   call PAR_MAX    (solve_sol(isolv) % cpu_dot(7))
                   call PAR_MAX    (solve_sol(isolv) % cpu_dot(9))
                   call PAR_MAX    (flops_dot_max)
                   !
                   call PAR_AVERAGE(6_ip,solve_sol(isolv) % cputi)
                   !
                   ! Units
                   !
                   call maths_unit(flops_max,unit_char,unit_factor)
                   flops_max = flops_max * unit_factor
                   flops_ave = flops_ave * unit_factor
                   call maths_unit(flops_dot_max,unit_dot_char,unit_dot_factor)
                   flops_dot_max = flops_dot_max * unit_dot_factor
                   flops_dot_ave = flops_dot_ave * unit_dot_factor
                   !
                   ! Write info
                   !
                   if( INOTSLAVE ) then

                      ii=len(trim(solve_sol(isolv) % wprob),kind=ip)+1_ip
                      write(nunit4,3700) &
                           trim(solve_sol(isolv) % wprob),&
                           woutl(1:ii),&
                           trim(solve_sol(isolv) % wsolv),&
                           solve_sol(isolv) % nsolv
                      if( solve_sol(isolv) % kfl_algso == 2 ) then
                         write(nunit4,3702)&
                              'NUMBER OF GROUPS=     ',&
                              solve_sol(isolv) % ngrou
                      end if

                      if( solve_sol(isolv) % kfl_algso >= 1 .and. solve_sol(isolv) % kfl_algso /= SOL_SOLVER_SPARSE_DIRECT ) then
                         write(nunit4,3703)&
                              'PRECONDITIONER:       ',&
                              trim(solve_sol(isolv) % wprec)
                         if( solve_sol(isolv) % kfl_preco == SOL_LINELET ) then
                            write(nunit4,3704)&
                                 'TOLERANCE PREC.=      ',&
                                 solve_sol(isolv) % toler
                         end if
                      end if

                      if( solve_sol(isolv) % kfl_algso >= 1 .and. solve_sol(isolv) % kfl_algso /= SOL_SOLVER_SPARSE_DIRECT ) then

                         write(nunit4,3701)&
                              solve_sol(isolv) % itsol(1),&
                              solve_sol(isolv) % itsol(2),&
                              solve_sol(isolv) % itsol(3)/solve_sol(isolv) % nsolv,&
                              solve_sol(isolv) % itsol(3),&
                              solve_sol(isolv) % cputi(1)

                         if( solve_sol(isolv) % kfl_schur > 0 ) then
                            !
                            ! Schur complement solver
                            !
                            routp( 1)  = sum(solve_sol(isolv) % cpu_schur)
                            routp( 1)  = max(zeror,routp(1))
                            dummr      = 100.0_rp / routp(1)
                            routp( 2)  = solve_sol(isolv) % cpu_schur(1)
                            routp( 3)  = solve_sol(isolv) % cpu_schur(1)*dummr       ! Invert Aii
                            routp( 4)  = solve_sol(isolv) % cpu_schur(2)
                            routp( 5)  = solve_sol(isolv) % cpu_schur(2)*dummr       ! Compute preconditioner P
                            routp( 6)  = solve_sol(isolv) % cpu_schur(3)
                            routp( 7)  = solve_sol(isolv) % cpu_schur(3)*dummr       ! Factorize preconditioner P
                            routp( 8)  = solve_sol(isolv) % cpu_schur(4)
                            routp( 9)  = solve_sol(isolv) % cpu_schur(4)*dummr       ! Solver initialization
                            routp(10)  = solve_sol(isolv) % cpu_schur(6)
                            routp(11)  = solve_sol(isolv) % cpu_schur(6)*dummr       ! Main loop: q = A p = Abb p - Abi Aii^{-1} Aib p
                            routp(12)  = solve_sol(isolv) % cpu_schur(7)
                            routp(13)  = solve_sol(isolv) % cpu_schur(7)*dummr       ! Main loop: x = x + alpha*p
                            routp(14)  = solve_sol(isolv) % cpu_schur(8)
                            routp(15)  = solve_sol(isolv) % cpu_schur(8)*dummr       ! Main loop: L z = q
                            routp(16)  = solve_sol(isolv) % cpu_schur(5)&
                                 -(solve_sol(isolv) % cpu_schur(6) &
                                 + solve_sol(isolv) % cpu_schur(7) &
                                 + solve_sol(isolv) % cpu_schur(8)  )
                            routp(16)  = max(zeror,routp(16))
                            routp(17)  = routp(16)*dummr                             ! Main loop: others
                            write(nunit4,3706) &
                                 routp( 1) ,  &
                                 routp( 2),routp( 3),routp( 4),routp( 5) , &
                                 routp( 6),routp( 7),routp( 8),routp( 9) , &
                                 routp(10),routp(11),routp(12),routp(13) , &
                                 routp(14),routp(15),routp(16),routp(17)
                         else
                            !
                            ! Iterative solver
                            !
                            routp     = 0.0_rp
                            routp(1)  = 100.0_rp / max(solve_sol(isolv) % cputi(1),zeror)
                            routp(2)  = solve_sol(isolv) % cputi(2) * routp(1)
                            routp(3)  = solve_sol(isolv) % cputi(3) * routp(1)
                            routp(4)  = solve_sol(isolv) % cputi(4) * routp(1)
                            routp(5)  = solve_sol(isolv) % cputi(5) * routp(1)
                            routp(6)  = solve_sol(isolv) % cputi(6) * routp(1)
                            dummr     = solve_sol(isolv) % cputi(1)-sum(solve_sol(isolv) % cputi(2:10))
                            dummr     = max(dummr,0.0_rp)
                            routp(10) = dummr * routp(1)
                            write(nunit4,3705)&
                                 solve_sol(isolv) % cputi(1),&
                                 solve_sol(isolv) % cputi(2),routp(2),&
                                 solve_sol(isolv) % cputi(3),routp(3),&
                                 solve_sol(isolv) % cputi(4),routp(4),&
                                 solve_sol(isolv) % cputi(5),routp(5)
                            if( solve_sol(isolv) % kfl_full_rows /= 0 ) then
                               write(nunit4,3708)&
                                    solve_sol(isolv) % cputi(6),routp(6)
                            end if
                            write(nunit4,3709)&
                                 dummr,routp(10)
                            if( solve_sol(isolv) % num_spmv > 0 ) then
                               write(nunit4,3707)&
                                    solve_sol(isolv) % num_spmv,&
                                    solve_sol(isolv) % cpu_spmv(4),&
                                    solve_sol(isolv) % cpu_spmv(5),&
                                    solve_sol(isolv) % cpu_spmv(4)/(solve_sol(isolv) % cpu_spmv(5)+zeror),&
                                    solve_sol(isolv) % cpu_spmv(6),&
                                    solve_sol(isolv) % cpu_spmv(7),&
                                    solve_sol(isolv) % cpu_spmv(6)/(solve_sol(isolv) % cpu_spmv(7)+zeror),&
                                    solve_sol(isolv) % cpu_spmv(8),&
                                    solve_sol(isolv) % cpu_spmv(9),&
                                    solve_sol(isolv) % cpu_spmv(8)/(solve_sol(isolv) % cpu_spmv(9)+zeror),&
                                    flops_ave,unit_char,&
                                    flops_max,unit_char
                            end if
                            if( solve_sol(isolv) % num_dot > 0 ) then
                               write(nunit4,3717)&
                                    solve_sol(isolv) % num_dot,&
                                    solve_sol(isolv) % cpu_dot(4),&
                                    solve_sol(isolv) % cpu_dot(5),&
                                    solve_sol(isolv) % cpu_dot(4)/(solve_sol(isolv) % cpu_dot(5)+zeror),&
                                    solve_sol(isolv) % cpu_dot(6),&
                                    solve_sol(isolv) % cpu_dot(7),&
                                    solve_sol(isolv) % cpu_dot(6)/(solve_sol(isolv) % cpu_dot(7)+zeror),&
                                    solve_sol(isolv) % cpu_dot(8),&
                                    solve_sol(isolv) % cpu_dot(9),&
                                    solve_sol(isolv) % cpu_dot(8)/(solve_sol(isolv) % cpu_dot(9)+zeror),&
                                    flops_dot_ave,unit_dot_char,&
                                    flops_dot_max,unit_dot_char
                            end if
                         end if
                      end if
                      write(nunit4,*)
                   end if
                end if
             end if
          end do

       case(38_ip)
          !
          ! Partition summary
          !
          !write(nunit4,'(a)') '* MPI partitioning'
          write(UNIT=nunit4,FMT=3800,IOSTAT=istat4) &
               (ioutp(ii),ii= 1,12),&
               routp(1),&
               (ioutp(ii),ii=13,17),&
               routp(2),&
               (ioutp(ii),ii=18,37)
          !if(istat4/=0) print*,'ERROR=',istat4

       case(39_ip)
          !
          ! Algebraic solver header
          !
          if( solve_sol(1) % kfl_solve == 1 ) then

             value_mix(1) = -solve_sol(1) % nzmat
             value_mix(2) =  solve_sol(1) % nzmat
             value_mix(3) = -solve_sol(1) % omp_chunk_size 
             value_mix(4) =  solve_sol(1) % omp_chunk_size 

             value_ave(1) = solve_sol(1) % nzmat
             value_ave(2) = solve_sol(1) % omp_chunk_size

             call PAR_MAX(4_ip,value_mix)
             call PAR_AVERAGE(1_ip,value_ave(1:1))
             call PAR_SUM(1_ip,value_ave(2:2))

             nzmat_min          = -value_mix(1)
             nzmat_max          =  value_mix(2)
             nzmat_ave          =  value_ave(1)
             omp_chunk_size_min = -value_mix(3) 
             omp_chunk_size_max =  value_mix(4) 
             omp_chunk_size_ave =  value_ave(2)     

             omp_scheduling = 0

             if(      solve_sol(1) % omp_schedule == SOL_OMP_OFF ) then
                omp_scheduling(1) = 1
             else if( solve_sol(1) % omp_schedule == SOL_OMP_STATIC ) then
                omp_scheduling(2) = 1
             else if( solve_sol(1) % omp_schedule == SOL_OMP_GUIDED ) then
                omp_scheduling(3) = 1
             else if( solve_sol(1) % omp_schedule == SOL_OMP_DYNAMIC ) then
                omp_scheduling(4) = 1
             end if
             call PAR_SUM(4_ip,omp_scheduling)
             if( omp_scheduling(4) > 0 ) omp_chunk_size_ave = omp_chunk_size_ave / omp_scheduling(4)

             if( INOTSLAVE ) then
                if( kfl_rstar == 2 ) then
                   write(solve_sol(1) % lun_solve,3929)
                end if
                if( solve_sol(1) % kfl_schur == 1 ) then
                   write(solve_sol(1) % lun_solve,3903) trim(solve_sol(1) % wsolv)
                   if( solve_sol(1) % kfl_scaii == 0 ) then
                      write(solve_sol(1) % lun_solve,3926) 'CHOLESKY IN SKYLINE FORMAT'
                   else
                      write(solve_sol(1) % lun_solve,3926) 'DIRECT SPARSE CSR'
                   end if
                else
                   write(solve_sol(1) % lun_solve,3903) trim(solve_sol(1) % wsolv)
                end if
                if( solve_sol(1) % kfl_blogs == 1 ) then
                   write(solve_sol(1) % lun_solve,3928) 'ON'
                end if
                if(      solve_sol(1) % kfl_algso ==  0 ) then
                else if( solve_sol(1) % kfl_algso == -1 ) then
                else if( solve_sol(1) % kfl_algso == -2 ) then
                else if( solve_sol(1) % kfl_algso == -3 ) then
                else

                   if(      solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
                        &   solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
                      !
                      ! Deflation
                      !
                      write(solve_sol(1) % lun_solve,3908) solve_sol(1) % ngrou
                      write(solve_sol(1) % lun_solve,3911) trim(direct_solver_name(solve_sol(1) % direct_solver_Deflation))

                   else if( solve_sol(1) % kfl_algso == SOL_SOLVER_GMRES ) then
                      !
                      ! GMRES
                      !
                      write(solve_sol(1) % lun_solve,3905) solve_sol(1) % nkryd
                      if( solve_sol(1) % kfl_ortho == 0 ) then
                         write(solve_sol(1) % lun_solve,3950) 'CLASSICAL GRAM-SCHMIDT'
                      else if( solve_sol(1) % kfl_ortho == 1 ) then
                         write(solve_sol(1) % lun_solve,3950) 'MODIFIED GRAM-SCHMIDT'
                      end if

                   end if

                   write(solve_sol(1) % lun_solve,3906) solve_sol(1) % solco
                   write(solve_sol(1) % lun_solve,3907) solve_sol(1) % miter
                   if( solve_sol(1) % kfl_adres == 0 ) then
                      write(solve_sol(1) % lun_solve,3918)
                   else
                      write(solve_sol(1) % lun_solve,3919) solve_sol(1) % adres
                   end if

                   memor_value   = int(nzmat_max,8) * 8_8
                   call memory_unit(memor_value,memor_char,memor_factor)
                   mem_nzmat_min = real(nzmat_min,rp) * memor_factor * 8_rp
                   mem_nzmat_max = real(nzmat_max,rp) * memor_factor * 8_rp
                   mem_nzmat_ave = real(nzmat_ave,rp) * memor_factor * 8_rp

                   if( ISEQUEN ) then
                      write(solve_sol(1) % lun_solve,3910) nzmat_max,rp,mem_nzmat_max,memor_char
                   else
                      write(solve_sol(1) % lun_solve,3915) nzmat_min,rp,mem_nzmat_min,memor_char
                      write(solve_sol(1) % lun_solve,3916) nzmat_max,rp,mem_nzmat_max,memor_char
                      write(solve_sol(1) % lun_solve,3917) nzmat_ave,rp,mem_nzmat_ave,memor_char
                   end if
                   if(     solve_sol(1) % kfl_format == SOL_CSR_FORMAT ) then
                      write(solve_sol(1) % lun_solve,3935) 'CSR'
                   else if( solve_sol(1) % kfl_format == SOL_COO_FORMAT ) then
                      write(solve_sol(1) % lun_solve,3935) 'COO'
                   else if( solve_sol(1) % kfl_format == SOL_ELL_FORMAT ) then
                      write(solve_sol(1) % lun_solve,3935) 'ELL'
                   end if
                   !
                   ! Parallelization
                   !
                   write(solve_sol(1) % lun_solve,3941) 
                   if( solve_sol(1) % kfl_full_rows == 0 ) then
                      write(solve_sol(1) % lun_solve,3942) 'PARTIAL ROW (SQUARE MATRIX)'
                   else
                      write(solve_sol(1) % lun_solve,3942) 'FULL ROW (RECTANGULAR MATRIX)'
                   end if
                   if( kfl_async == 0 ) then
                      write(solve_sol(1) % lun_solve,3936) 'BLOCKING'
                   else
                      write(solve_sol(1) % lun_solve,3936) 'NON-BLOCKING'
                   end if
                   if( par_hybrid /= 0 ) then
                      if( omp_scheduling(1) > 0 ) then
                         write(solve_sol(1) % lun_solve,3943) 'OFF',omp_scheduling(1)
                      end if
                      if( omp_scheduling(2) > 0 ) then
                         write(solve_sol(1) % lun_solve,3943) 'STATIC SCHEDULING',omp_scheduling(2)
                      end if
                      if( omp_scheduling(3) > 0 ) then
                         write(solve_sol(1) % lun_solve,3943) 'GUIDED SCHEDULING',omp_scheduling(3)
                      end if
                      if( omp_scheduling(4) > 0) then
                         write(solve_sol(1) % lun_solve,3943) 'DYNAMIC SCHEDULING',omp_scheduling(4)
                         write(solve_sol(1) % lun_solve,3944) omp_chunk_size_min                 
                         write(solve_sol(1) % lun_solve,3945) omp_chunk_size_max                 
                         write(solve_sol(1) % lun_solve,3946) omp_chunk_size_ave                 
                      end if
                   end if
                end if
                !
                ! Coarse solver
                !
                if( solve_sol(1) % kfl_coarse /= 0 ) then
                   write(solve_sol(1) % lun_solve,3913)
                   write(solve_sol(1) % lun_solve,3914) solve_sol(1) % ngrou
                   write(solve_sol(1) % lun_solve,3932) trim(direct_solver_name(solve_sol(1) % direct_solver_coarse))
                end if
                !
                ! Preconditioner header
                !
                if(  solve_sol(1) % kfl_algso /= 14 .and. &
                     solve_sol(1) % kfl_algso /= 20 .and. &
                     solve_sol(1) % kfl_algso /= 21 .and. &
                     solve_sol(1) % kfl_algso /= 22 ) then
                   write(solve_sol(1) % lun_solve,3920) solve_sol(1) % wprec
                   if(    solve_sol(1) % kfl_algso == SOL_SOLVER_CG           .or. &
                        & solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG  .or. &
                        & solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_CG ) then
                      write(solve_sol(1) % lun_solve,3934) 'M^{-1/2} A M^{-1/2} M^{1/2} x = M^{-1/2} b'
                   else
                      if( solve_sol(1) % kfl_leftr == 0 ) then
                         write(solve_sol(1) % lun_solve,3931) 'M^{-1} A x = M^{-1} b'
                      else
                         write(solve_sol(1) % lun_solve,3933) 'A M^{-1} M x = b'
                      end if
                   end if
                   if( solve_sol(1) % kfl_preco == 4 ) then
                      write(solve_sol(1) % lun_solve,3909) &
                           solve_sol(1) % toler,solve_sol(1) % nline,solve_sol(1) % nlin1,solve_sol(1) % npoin
                   end if
                   if( solve_sol(1) % kfl_preco == SOL_MULTIGRID ) then
                      write(solve_sol(1) % lun_solve,3922) solve_sol(1) % ngrou
                      if( solve_sol(1) % kfl_defso == 0 ) then
                         if( solve_sol(1) % kfl_defas == 0 ) then
                            write(solve_sol(1) % lun_solve,3921) 'CHOLESKY (SKYLINE)'
                         else
                            write(solve_sol(1) % lun_solve,3921) 'DIRECT SPARSE (CSR)'
                         end if
                      else
                         write(solve_sol(1) % lun_solve,3921) 'ITERATIVE SPARSE (CSR)'
                      end if
                      if( solve_sol(1) % kfl_defpr == 2 ) then
                         write(solve_sol(1) % lun_solve,3924) 'JACOBI WITH DIAGONAL PRECONDITIONER'
                      else if( solve_sol(1) % kfl_defpr == 4 ) then
                         write(solve_sol(1) % lun_solve,3923) 'JACOBI WITH LINELET PRECONDITIONER'
                         write(solve_sol(1) % lun_solve,3909) &
                              solve_sol(1) % toler,solve_sol(1) % nline,solve_sol(1) % nlin1,solve_sol(1) % npoin
                      end if
                      write(solve_sol(1) % lun_solve,3925) solve_sol(1) % itpre
                   end if
                   if( solve_sol(1) % kfl_preco == SOL_ABB ) then
                      if( solve_sol(1) % kfl_scpre == 0 ) then
                         write(solve_sol(1) % lun_solve,3927) 'CHOLESKY IN SKYLINE FORMAT'
                      else
                         write(solve_sol(1) % lun_solve,3927) 'DIRECT SPARSE CSR'
                      end if
                   end if
                   !
                   ! RAS
                   !
                   if( solve_sol(1) % kfl_preco == SOL_RAS ) then
                      if( solve_sol(1) % kfl_block_ras == 0 ) then
                         write(solve_sol(1) % lun_solve,3932) trim(direct_solver_name(solve_sol(1) % direct_solver_RAS(1)))
                      else
                         write(solve_sol(1) % lun_solve,3932) trim(direct_solver_name(solve_sol(1) % direct_solver_RAS(1)))//' (BY BLOCK)'
                      end if
                   end if
                   !
                   ! RAS
                   !
                   if( solve_sol(1) % kfl_preco == SOL_AII ) then
                      write(solve_sol(1) % lun_solve,3932) trim(direct_solver_name(solve_sol(1) % direct_solver_Block_LU))
                   end if

                end if
                !
                ! Solver header
                !
                if( solve_sol(1) % kfl_algso/=0.and.solve_sol(1) % kfl_algso/=-1)&
                     write(solve_sol(1) % lun_solve,3930)
                if( solve_sol(1) % kfl_cvgso == 1 ) then
                   write(solve_sol(1) % lun_cvgso,3940)
                end if
             end if
          end if

       case(40_ip)
          !
          ! User and system info
          !
          if( present(CHA_LIST) ) then
             write(nunit4,4100)
             write(nunit4,4102)
             do ii = 1,size(CHA_LIST)
                write(nunit4,4101) adjustl(trim(CHA_LIST(ii)))
             end do
          end if

       case(41_ip)
          !
          ! Makefile info
          !
          if( present(CHA_LIST) ) then
             write(nunit4,4103)
             do ii = 1,size(CHA_LIST)
                write(nunit4,4101) adjustl(trim(CHA_LIST(ii)))
             end do
          end if

       case(42_ip)
          !
          ! Boundary condition codes
          !
          write(nunit4,4200)

       case(43_ip)
          !
          ! List of boundary condition codes
          !
          !write(nunit4,4300) (ioutp(ii),ii=1,ioutp(49))
          !write(nunit4,4301) ioutp(50)
          if( ioutp(2) == -(mcodb+1) ) then
             write(nunit4,4301) ioutp(1)
          else if( ioutp(3) == -(mcodb+1) ) then
             write(nunit4,4302) ioutp(1),ioutp(2)
          else
             write(nunit4,4303) ioutp(1),ioutp(2),ioutp(3)
          end if

       case(44_ip)
          !
          ! Geomtrical condition nodes
          !
          write(nunit4,4400) (ioutp(ii),ii=1,10)

       case(45_ip)
          !
          ! Bad edge
          !
          write(nunit4,100) 'BAD EDGES HAVE BEEN FOUND'
          write(nunit4,4500) ioutp(1),ioutp(2),ioutp(3),ioutp(4)
          call runend('BAD EDGES HAVE FOUND: SEE LOG FILE')

       case(46_ip)
          !
          ! Sets results: time header
          !
          write(nunit4,4600) ioutp(1),routp(1)

       case(47_ip)
          !
          ! Error in kernel
          !
          if( present(messa) ) call runend(trim(messa))

       case(48_ip)
          !
          ! Error in kernel
          !
          if( present(messa) ) then
             if( live_info_config%lun_livei /= 0 ) call livinf(0_ip,trim(messa),0_ip)
             call runend(trim(messa))
          end if

       case(49_ip)
          !
          ! CPU times (subelm)
          !
          write(nunit4,4900) &
                                !routp( 1) ,  &
               routp( 2),routp( 3),routp( 4),routp( 5) , &
               routp( 6),routp( 7),routp( 8),routp( 9) , &
               routp(10),routp(11),routp(12),routp(13) , &
               routp(14),routp(15),routp(16),routp(17) , &
               routp(18),routp(19)

       case(50_ip)
          !
          ! Witness points
          !
          !write(nunit4,'(a)') '* Witness points'
          write(nunit4,5000)

       case(51_ip)
          !
          ! Witness points: List of missing witness points
          !
          if( ioutp(2) == 0 ) write(nunit4,5101)
          write(nunit4,5100) ioutp(1)

       case(52_ip)
          !
          ! Witness points: OK
          !
          write(nunit4,5200)

       case(53_ip)
          !
          ! Lagrangian particles: results
          !
          !write(nunit,5300)
          !write(nunit,5301)
          write(nunit4,5302)


       case(54_ip)
          !
          ! Lagrangian particles: convergence
          !
          write(nunit4,5400)

       case(55_ip)
          !
          ! Lagrangian particles: convergence with thermodynamic data
          !
          write(nunit4,5500)

       case(56_ip)
          !
          ! Header CPU times of other operations
          !
          write(nunit4,5600)

       case(57_ip)
          !
          ! Lagrangian particles: deposition map
          !
          write(nunit4,5700)
          !write(nunit4,5701)

       case(58_ip)
          !
          ! Support geometry projection
          !
          write(nunit4,5800) routp(1),routp(2),routp(3)

       case(59_ip)
          !
          ! Wake model
          !
          write(nunit4,5900)

       case(60_ip)
          !
          ! Wake model
          !
          write(nunit4,6000) ioutp(1),routp(1:13)

       case(61_ip)
          !
          ! Generic message
          !
          write(nunit4,6100) coutp(1),ioutp(1)

       case(62_ip)
          !
          ! Generic message
          !
          write(nunit4,6200) coutp(1),routp(1)

       case(63_ip)
          !
          ! Header zone-wise partition
          !
          write(nunit4,6300)

       case(64_ip)
          !
          ! Zone-wise partition
          !
          write(nunit4,6400) ioutp(1),ioutp(2)

       case(65_ip)
          !
          ! Thrust and power
          !
          write(abs(nunit4),6500) ioutp(1),ioutp(2), routp(1),routp(2),routp(3), routp(4), routp(5), routp(6), routp(7)

       case(66_ip)
          !
          ! Inter-color Communicators
          !
          !write(nunit4,'(a)') '* Communicators'
          write(nunit4,6600)
          write(nunit4,6601) ioutp(1),ioutp(2),ioutp(3),&
               &            ioutp(4),ioutp(5)

       case(67_ip)
          !
          ! List of inter-color Communicators
          !
          coutp(1) = 'CODE= '
          coutp(2) = 'ZONE= '
          coutp(3) = 'SUBD= '
          !write(nunit4,6700) trim(coutp(1)),ioutp(1),&
          !     &            trim(coutp(2)),ioutp(2),&
          !     &            trim(coutp(3)),ioutp(3),&
          !     &            trim(coutp(1)),ioutp(4),&
          !     &            trim(coutp(2)),ioutp(5),&
          !     &            trim(coutp(3)),ioutp(6)
          write(nunit4,6700) trim(coutp(1)),ioutp(1),ioutp(4),&
               &            trim(coutp(2)),ioutp(2),ioutp(5),&
               &            trim(coutp(3)),ioutp(3),ioutp(6)

       case(68_ip)
          !
          ! Parall bin structure
          !
          !write(nunit4,'(a)') '* Parallel bin structure'
          write(nunit4,6800) ioutp(1),ioutp(2),ioutp(3)

       case(69_ip)
          !
          ! Max slave memory (outmem)
          !
          write(nunit4,6900) &
               routp(1),trim(coutp(1)),&
               routp(2),trim(coutp(2))

       case(70_ip)
          !
          ! Load balance and communications
          !
          write(nunit4,7000)

       case(71_ip)
          !
          ! Load balance and communications
          !
          write(nunit4,7100) routp(1),routp(2),routp(3)

       case(72_ip)
          !
          ! Hybrid parallelization
          !
          if( INOTSLAVE ) write(nunit4,7200)
          dumm1 = par_hybrid
          dumm2 = par_hybrid
          call PAR_MAX(dumm1)
          call PAR_MIN(dumm2)
          if( dumm1 == 0 ) then
             if( INOTSLAVE ) write(nunit4,7206)
          else
             ioutp(1) = par_omp_num_threads
             ioutp(2) = par_omp_num_blocks
             ioutp(3) = par_omp_nelem_chunk
             ioutp(4) = par_omp_nboun_chunk
             ioutp(5) = par_omp_npoin_chunk
             if( INOTMASTER ) then
                if( par_hybrid == PAR_OPENMP_COLORING ) then
                   ioutp(6) = par_omp_num_colors
                else if( par_hybrid == PAR_OMPSS ) then
                   ioutp(6) = size(ompss_domains,kind=ip)
                   ioutp(7) = 0
                   do ii = 1,num_subd_par        
                      ioutp(7) = ioutp(7) + size(ompss_domains(ii) % neighbours,kind=ip)
                   end do
                   ioutp(7) = ioutp(7) / max(1_ip,num_subd_par)
                else
                   ioutp(6) = 0
                   ioutp(7) = 0
                end if
             end if

             call PAR_AVERAGE(5_ip,ioutp(3:7),'IN MY CODE')

             if( INOTSLAVE ) then
                if( dumm1 /= dumm2 ) then
                   coutp(1) = 'CO-EXECUTION WITH DIFFERENT STRATEGIES... NUMBERS DONT MAKE SENSE'
                else if( par_hybrid == PAR_OPENMP_NO_COLORING ) then
                   coutp(1) = 'OPENMP WITHOUT COLORING'
                else if( par_hybrid == PAR_OPENMP_COLORING ) then
                   coutp(1) = 'OPENMP WITH COLORING'
                   if( par_omp_coloring_alg == 0 ) then
                      coutp(1) = trim(coutp(1)) // ' (GREEDY ALGORITHM)'
                   else if( par_omp_coloring_alg == 1 ) then
                      coutp(1) = trim(coutp(1)) // ' (HOMEMADE ALGORITHM)'
                   end if
                else if( par_hybrid == PAR_OMPSS ) then
                   coutp(1) = 'OMPSS WITH MULTI-DEPENDENCIES'
                end if
                coutp(2) = 'OPENMP WITHOUT COLORING'
#ifdef ALYA_DLB
                coutp(3) = 'ON'
#else
                coutp(3) = 'OFF'
#endif
                write(nunit4,7205) &
                     trim(coutp(1)),trim(coutp(2)),trim(coutp(3)),&
                     ioutp(1),ioutp(2),par_omp_granularity,ioutp(3:5)
                if( par_hybrid == PAR_OPENMP_COLORING .and. par_omp_num_colors > 0 ) then
                   write(nunit4,7202) ioutp(6)
                   !write(nunit4,7203)
                   !do ii = 1,par_omp_num_colors
                   !   write(nunit4,7201) ii,par_omp_ia_colors(ii+1)-par_omp_ia_colors(ii)
                   !end do
                else if( par_hybrid == PAR_OMPSS ) then
                   write(nunit4,7204) ioutp(6),ioutp(7)
                end if
             end if
          end if
          !
          ! Co-execution
          !
          dummi = 0
#ifdef OPENACCHHH
          dummi = 1
#endif
          call PAR_SUM(dummi)      
          if( dummi > 0 .and. INOTSLAVE ) then
             write(nunit4,7210)
             write(nunit4,7211) npart-dummi,dummi
          end if

       case(73_ip)
          !
          ! Optimization options
          !
          !write(lun_outpu4,'(a)') '* Optimization'
#ifdef MPI3
          wmpi3 = 'ON'
#else
          wmpi3 = 'OFF'
#endif
#ifdef BLAS
          wblas = 'ON'
#else
          wblas = 'OFF'
#endif
          if( IMASTER ) write(lun_outpu4,7300) wmpi3,wblas
          if( IMASTER ) then
          write(lun_outpu4,7307) VECTOR_SIZE        
          end if

          if( kfl_element_to_csr == 0 ) then
             if( IMASTER ) write(lun_outpu4,7301)
          else
             memor_value4 = size(lezdo,kind=ip)
             call PAR_AVERAGE(memor_value4)
             memor_value = int(ip,8) * int(memor_value4,8)
             call memory_unit(memor_value,memor_char,memor_factor)
             if( IMASTER ) write(lun_outpu4,7302) real(memor_value,rp) * memor_factor,memor_char
          end if

          if( kfl_savda == 0 ) then
             if( IMASTER ) write(lun_outpu4,7303)
          else
             call element_data_base_memory(memor_value4)
             call PAR_AVERAGE(memor_value4)
             memor_value = int(memor_value4,8)
             call memory_unit(memor_value,memor_char,memor_factor)
             if( IMASTER ) write(lun_outpu4,7304) real(memor_value,rp) * memor_factor,memor_char
          end if

#ifdef OPENACCHHH
          ngpu = 1
#else
          ngpu = 0
#endif
          call PAR_SUM(ngpu)
          if( IMASTER ) then
             if( ngpu > 0 ) then
                write(lun_outpu4,7308) ngpu,npart-ngpu
             end if

             if(      kfl_renumbering_npoin == 0 ) then
                me100 = 'NONE'
             else if( kfl_renumbering_npoin == 1 ) then
#if   defined V5METIS 
                me100 = 'METIS 5.0.2'
#elif defined V51METIS 
                me100 = 'METIS 5.1.0'
#else
                me100 = 'METIS 4.0'
#endif
             else if( kfl_renumbering_npoin == 2 ) then
                me100 = 'SFC'
             else if( kfl_renumbering_npoin == 3 ) then
                me100 = 'CUTHILL-MCKEE'
             endif
             write(lun_outpu4,7309) trim(me100)

             if( kfl_renumbering_nelem == 0 ) then
                me100 = 'NONE'
             else if( kfl_renumbering_npoin == 1 ) then
                me100 = 'MESH GRAPH FIRST TOUCH'
             end if
             write(lun_outpu4,7310) trim(me100)
          end if

       case(74_ip)
          !
          ! Header direct solvers
          !
          ipass = 0
          !write(lun_outpu4,7400)

       case(75_ip)
          !
          ! Write direct solvers info
          !
          if( ipass == 0 ) then
             ipass = 1
             write(lun_outpu4,7400)
          end if
          write(lun_outpu4,7500)&
               trim(coutp(1)),trim(coutp(2)),trim(coutp(3)),trim(coutp(4)),&
               ioutp(4),&
               ioutp(1),routp(1),routp(4),routp(7),&                     ! Initialization
               ioutp(2),routp(2),routp(5),routp(8),routp(12),routp(13),& ! Factorization
               ioutp(3),routp(3),routp(6),routp(9),&                     ! Solution
               routp(10),trim(coutp(5)),&
               routp(11),trim(coutp(5))

       case(76_ip)
          !
          ! Typical operations
          !
          !write(lun_outpu4,'(a)') '* Timings algebraic solver operations'
          !
          ! Typical operations
          !
          !write(lun_outpu4,'(a)') '* Timings algebraic solver operations'
          write(lun_outpu4,7600) & 
               kfl_algebra_operations, &
               routp( 1),trim(coutp(1)),routp( 2),trim(coutp(1)),routp( 3),&
               routp( 4),trim(coutp(1)),routp( 5),trim(coutp(1)),routp( 6),&
               routp( 7),trim(coutp(1)),routp( 8),trim(coutp(1)),routp( 9)

          if( par_omp_num_threads /= 0 ) &
               write(lun_outpu4,7601)  &
               routp(10),trim(coutp(1)),routp(11),trim(coutp(1)),routp(12),&
               routp(13),trim(coutp(1)),routp(14),trim(coutp(1)),routp(15),&
               routp(16),trim(coutp(1)),routp(17),trim(coutp(1)),routp(18)

#ifdef NINJA
          write(lun_outpu4,7604)  &
               routp(37),trim(coutp(3)),routp(38),trim(coutp(3)),routp(39),&
               routp(40),trim(coutp(3)),routp(41),trim(coutp(3)),routp(42),&
               routp(43),trim(coutp(3)),routp(44),trim(coutp(3)),routp(45)        
#endif


          write(lun_outpu4,7602)  &
               routp(19),trim(coutp(2)),routp(20),trim(coutp(2)),routp(21),&
               routp(22),trim(coutp(2)),routp(23),trim(coutp(2)),routp(24),&
               routp(25),trim(coutp(2)),routp(26),trim(coutp(2)),routp(27)

          if( par_omp_num_threads /= 0 ) &
               write(lun_outpu4,7603)  &
               routp(28),trim(coutp(2)),routp(29),trim(coutp(2)),routp(30),&
               routp(31),trim(coutp(2)),routp(32),trim(coutp(2)),routp(33),&
               routp(34),trim(coutp(2)),routp(35),trim(coutp(2)),routp(36)

#ifdef NINJA
          write(lun_outpu4,7605)  &
               routp(46),trim(coutp(4)),routp(47),trim(coutp(4)),routp(48),&
               routp(49),trim(coutp(4)),routp(50),trim(coutp(4)),routp(51),&
               routp(52),trim(coutp(4)),routp(53),trim(coutp(4)),routp(54)        
#endif

       case(77_ip)
          !
          ! Partitioning
          !
          if( .not. READ_AND_RUN() ) then
             if(      kfl_partition_par == PAR_METIS4 ) then
#if   defined V5METIS 
                me100 = 'METIS 5.0.2'
#elif defined V51METIS 
                me100 = 'METIS 5.1.0'
#else
                me100 = 'METIS 4.0'
#endif
                write(lun_outpu4,7700) trim(me100)
                if( kfl_parti_par == 1 ) then
                   write(lun_outpu4,7705) 'NODES'
                else if( kfl_parti_par == 2 ) then
                   write(lun_outpu4,7705) 'FACES'
                end if
             else if( kfl_partition_par < 0 ) then
                write(lun_outpu4,7700) 'READ FROM ELEMENT FIELD '//trim(integer_to_string(-kfl_partition_par))
             else if( kfl_partition_par == PAR_SFC ) then
                write(lun_outpu4,7700) 'HILBERT SPACE FILLING CURVE'
             else if( kfl_partition_par == PAR_ZOLTAN ) then
                write(lun_outpu4,7700) 'ZOLTAN LIBRARY HSFC'
             else if( kfl_partition_par == PAR_ORIENTED_BIN ) then
                write(lun_outpu4,7700) 'ORIENTED BIN'
             end if
             !
             ! Sequential/parallel mode
             !
             if( kfl_parseq_par == PAR_SEQUENTIAL_PARTITION ) then
                write(lun_outpu4,7708) 'SEQUENTIAL'
             else
                write(lun_outpu4,7708) 'PARALLEL'
             end if
             !
             ! Weights
             !
             if( kfl_weigh_par == 0 ) then
                write(lun_outpu4,7706) 'GAUSS POINTS'
             else if( kfl_weigh_par < 0 ) then
                write(lun_outpu4,7706) 'NONE'
             else if( kfl_weigh_par > 0 ) then
                write(lun_outpu4,7706) 'FROM FIELD'
             end if
             !
             ! Timings
             !
             if( kfl_parseq_par == PAR_PARALLEL_PARTITION ) then
                dummr = sum(routp(1:5)+zeror)
                write(lun_outpu4,7709) &
                     dummr,&
                     (routp(ii),routp(ii)/dummr*100.0_rp,ii=1,5)
             end if
          end if

       case(78_ip)
          !
          ! Integration rule test header
          !
          write(lun_outpu4,7800)

       case(79_ip)
          !
          ! Integration rule test
          !
          write(lun_outpu4,7900) trim(coutp(1)),ioutp(1),ioutp(2),routp(1),routp(2),trim(coutp(2))

       case(80_ip)
          !
          ! Integration rule test header
          !
          write(lun_outpu4,8000)

       case(81_ip)
          !
          ! Integration rule test
          !
          if( ioutp(2) == 0 ) then
             write(lun_outpu4,8103) trim(coutp(1))
          else
             write(lun_outpu4,8102) trim(coutp(1)),trim(coutp(2)),ioutp(1),ioutp(2)
          end if

       case(82_ip)
          !
          ! Halos
          !
          write(lun_outpu4,8200) routp(1:6)

       case(83_ip)
          !
          ! Halos
          !
          dummr = 100.0_rp / (sum(routp(1:9))+zeror)
          write(nunit4,8300) &
               ioutp(1),&
               sum(routp(1:9)),&
               routp(1),routp(1)*dummr,&
               routp(2),routp(2)*dummr,&
               routp(3),routp(3)*dummr,&
               routp(4),routp(4)*dummr,&
               routp(5),routp(5)*dummr,&
               routp(6),routp(6)*dummr,&
               routp(7),routp(7)*dummr,&
               routp(8),routp(8)*dummr,&
               routp(9),routp(9)*dummr

       case(84_ip)
          !
          ! Boundary condition codes used
          !
          write(nunit4,8400) trim(coutp(1))

       case(85_ip)
          !
          ! Boundary condition codes used
          !
          write(nunit4,8500) trim(coutp(1))

       case(86_ip)
          !
          ! Boundary condition codes used
          !
          write(nunit4,8600) ioutp(1:8)

       case(87_ip)
          !
          ! Auto-tuning header
          !
          write(nunit4,8700) 

       case(88_ip)
          !
          ! Auto-tuning problem
          !
          write(nunit4,8800) trim(coutp(1)) 
       if( nunit4 /= 0_4 .and. INOTSLAVE ) call iofile_flush_unit(nunit4)
       case(89_ip)
          !
          ! Auto-tuning problem
          !
          if( ioutp(10) > 0 ) write(nunit4,8902) ioutp(10)   ! Guide
          if( ioutp(11) > 0 ) write(nunit4,8900) ioutp(11)   ! Static
          if( ioutp(12) > 0 ) write(nunit4,8904) ioutp(12)   ! Off
          if( ioutp(13) > 0 ) write(nunit4,8901) ioutp(13),& ! Dynamic
               ioutp(3),ioutp(4),ioutp(5)

          write(nunit4,8903) routp(3),trim(coutp(1)),&
               &            routp(4),trim(coutp(1)),&
               &            routp(5),trim(coutp(1)),&
               &            routp(6)

       case(90_ip)
          !
          ! Summary of computing time: title
          !
          write(nunit4,9000) trim(coutp(1)),max(routp(1),0.0_rp)

       case(91_ip)

          if(run_config%log) write(nunit4,9100)

       case(92_ip)

          write(nunit4,9200) ioutp(1)

       case(93_ip)

          write(nunit4,9300) trim(coutp(1))

       case(94_ip)

          write(nunit4,9400,advance='no') 
          write(nunit4,9401,advance='no') ioutp(1)

       case(95_ip)
          !
          ! Partitioning statistics
          !
          if( kfl_partition_par == PAR_SFC ) then

             write(lun_outpu4,7702) boxes_fine_par(1:ndime)
             write(lun_outpu4,7703) boxes_coarse_par(1:ndime)
             write(nunit4,9500) routp(5),routp(1:4)

          else if( kfl_partition_par == PAR_ORIENTED_BIN ) then
             write(lun_outpu4,7701) boxes_fine_par(1:ndime)
             if( ndime == 2 ) then
                write(lun_outpu4,7707) vect_partition_par(1:ndime)
             else
                write(lun_outpu4,7704) vect_partition_par(1:ndime)
             end if
          end if

       case(96_ip)
          !
          ! Element weights
          !
          write(nunit4,9600)

       case(97_ip)
          !
          ! Element weights
          !
          if( associated(INT_LIST) ) then
             write(nunit4,9700) namod(ioutp(1)),routp(1)
             do ii = 1,ioutp(2)                
                write(nunit4,9701) element_type(INT_LIST(ii)) % name,max(routp(ii+1),0.0_rp)
             end do
          end if

       case(98_ip)
          !
          ! Timings
          !
          write(nunit4,9800)

       case(99_ip)
          !
          ! Element assembly
          !
          call maths_unit(routp(10),unit_char,unit_factor)
          routp(10) = routp(10) * unit_factor ! Scale timings per element
          routp( 9) = min(routp(9),9999.0_rp) ! Limit variability to avoid overflow
          write(nunit4,9900)&
               trim(coutp(1)),&               
               routp( 1) , routp( 2)    ,   & ! Matrix onstruction average
               routp( 3) , routp( 4)    ,   & ! Matrix construction max
               routp( 5)                ,   & ! Load balance              
               routp(10),trim(unit_char),   & ! Average element cost per assembly
               routp( 6)                ,   & ! Theoretical Matrix onstruction average
               routp( 7)                ,   & ! Theoretical Matrix construction max
               routp( 8)                ,   & ! Theoretical Load balance
               routp( 9)                      ! Maximum variability %

       case(100_ip)
          !
          ! Ather assemblies
          !
          if( ioutp(1) == 1 ) then
             call maths_unit(routp(10),unit_char,unit_factor)
             routp(10) = routp(10) * unit_factor            
             write(nunit4,10000)&
                  trim(coutp(1)),&              
                  routp( 1) , routp( 2)     , & ! Matrix onstruction average
                  routp( 3) , routp( 4)     , & ! Matrix construction max
                  routp( 5) ,                 & ! Load balance
                  routp(10) , trim(unit_char)   ! Average cost per entity per assembly
          else if( ioutp(1) == 2 ) then
             dummi = 26-len_trim(coutp(1))
             if (dummi<0) call runend("mod_outfor: length of coutp(1) is > 26, coutp(1)='"//trim(coutp(1))//"'")
             write(nunit4,10002)&
                  trim(coutp(1))//':'//repeat(' ',dummi),&              
                  routp( 1) , routp( 2)         ! Max timing
          else if( ioutp(1) == 0 ) then
             write(nunit4,10001)&
                  trim(coutp(1)),&              
                  routp( 1) , routp( 2)     , & ! Matrix onstruction average
                  routp( 3) , routp( 4)     , & ! Matrix construction max
                  routp( 5)                     ! Load balance
          end if


       case(101_ip)
          !
          ! Deposition
          !
          write(nunit4,10100)

       case(102_ip)
          !
          ! Modules info
          !
          if( adjustl(coutp(1)(1:3))/='END' ) then
             write(nunit4,10200)
             ii=1
             do while(coutp(ii)/='END'.and.ii<10)
                write(nunit4,10201) adjustl(trim(coutp(ii)))
                ii=ii+1
             end do
          end if

       case(103_ip)
          !
          ! Repartitioning
          !
          if( ioutp(1) > 0 ) then
             routp(12) = routp(2) * 100.0_rp / (routp(1) + zeror)
             routp(13) = routp(3) * 100.0_rp / (routp(1) + zeror)
             routp(14) = routp(4) * 100.0_rp / (routp(1) + zeror)
             routp(15) = routp(5) * 100.0_rp / (routp(1) + zeror)
             write(nunit4,10300) & 
                  ioutp(1),&
                  routp(1),&
                  routp(2),routp(12),&
                  routp(3),routp(13),&
                  routp(4),routp(14),&
                  routp(5),routp(15)
          end if
          
       case(104_ip)
          !
          ! Parallel Performance title
          !
          if( routp(1) >= 0.0_rp ) then
             write(nunit4,10400) 
             write(nunit4,10401) adjustl('GLOBAL:       '),routp(1:7) 
          end if

       case(105_ip)
          !
          ! Parallel Performance of modules
          !
          if( routp(1) >= 0.0_rp ) then
             write(nunit4,10401) adjustl(trim(namod(ioutp(1))))//' MODULE:',routp(1:7)
          end if

       case(106_ip)
          !
          ! Material data header
          !
          write(nunit4,10600)

       case(107_ip)
          !
          ! Material data
          !
          write(nunit4,10700) ioutp(1),routp(1)

       case(108_ip)
          !
          ! Properties
          !      
          write(nunit4,10800) routp(1),routp(2)
             
       case(109_ip)
          !
          ! Details
          !      
          write(nunit4,10900)
          
       case(110_ip)
          !
          ! Details
          !      
          write(nunit4,11000) routp(1)
          
       case(111_ip)
          !
          ! Details
          !
          dummi = 26-len_trim(coutp(1))
          if (dummi<0) call runend("mod_outfor: length of coutp(1) is > 26, coutp(1)='"//trim(coutp(1))//"'")

          write(nunit4,10002)&
               trim(coutp(1))//':'//repeat(' ',dummi),&
               routp(1),&
               routp(2)    
          
       end select

       if( nunit4 /= 0_4 .and. INOTSLAVE ) call iofile_flush_unit(nunit4)

    end if
    !
    ! When slaves should stop
    !
    if( ISLAVE ) then
       if( itask == 4 ) then
          if( present(messa) ) then
             if( trim(messa) == '1' ) then
                call runend(trim(messa)//' ERROR HAS BEEN FOUND IN '//&
                     trim(namod(modul))//' DATA FILE')
             else
                call runend(trim(messa)//' ERRORS HAVE BEEN FOUND IN '//&
                     trim(namod(modul))//' DATA FILE')
             end if
          else
             call runend('ERRORS HAVE BEEN FOUND IN '//&
                  trim(namod(modul))//' DATA FILE')
          end if
       end if
    end if
    !
    ! Formats
    !
1   format(//,&
         & 5x,//,&
         & 5x,'--------------------------------------------------------',/,&
         & 5x,/,&
         & 5x,'|- AN ERROR HAS BEEN DETECTED:',/,/,&
         & 9x,a,/)
100 format(&
         & 5x,'|- ERROR:   ',a)
200 format(&
         & 5x,'|- WARNING: ',a)
201 format(&
         & 5x,'|- ! WARNING: ',a)
300 format(&
         & 5x,'|- ',a)
500 format(&
         & '# Time step: ',i8,', Outer iter:  ',i5,', Inner iter:  ',i5)
600 format(//,&
         & 5x,'|- END OF ',a,' RUN',/)
800 format(5x,&
         & 5x,' ',a6,':',/,&
         & 5x,'      - Critical time step dtc: ',es10.3,/,&
         & 5x,'      - Ratio dt/dtc:           ',es10.3 )
801 format(5x,&
         & 5x,' ',a6,':',/,&
         & 5x,'      - Critical time step dtc: ',e12.6,/,&
         & 5x,'      - Using local time step. min= ',e12.6,' max= ',e12.6)
900 format(&
         & 5x,'      Residual of ',a17,' : ',es10.3)
1000 format(/,&
         & 5x,'      Global iteration number: ',i5,', block number: ',i2,/,&
         & 5x,'      Current CPU time:        ',es10.3,/,&
         & 5x,'      Elapsed CPU time:        ',es10.3)
1100 format(///,&
         & 5x,'|- THIS IS A TIME RESTART RUN: CONTINUING... ')
1200 format(///,&
         & 5x,'   --------------------------',/,&
         & 5x,'|- ',a6,' OUTPUT FOR PROBLEM: ',a,/,&
         & 5x,'   --------------------------',/)
1400 format(//,&
         & 5x,'|- END OF ANALYSIS',/)
1500 format(/,&
         &    5x,'   |- TIME STEP NUMBER: ',i12,/&
         &    5x,'      Current time step dt = ',es10.3,/&
         &    5x,'      Current time t       = ',es10.3)
1700 format(//,&
         & 5x,'   -----------------------',/&
         & '*',4x,'|- SUMMARY OF PROBLEM DATA',/,&
         & 5x,'   -----------------------',/,&
         & /,&
         & /,&
         & 5x,'   |- GENERAL DATA:',/,&
         & /,&
         & 5x,'      Modules solved:       ',a,/,&
         & 5x,'      Blok definition:      ',a,/,&
         & 5x,'      Time step strategy:   ',a,/,&
         & 5x,'      Time step=            ',es11.4,/,&
         & 5x,'      Time range=           ','[',es11.4,',',es11.4,']',/,&
         & /,&
         & 5x,'   |- DOMAIN DATA:',/,&
         & /,&
         & 5x,'      Domain dimension=     ',i11,     /,&
         & 5x,'      --------------------- ',        /,&
         & 5x,'      Bandwith (average)=  ',i11,      /,&
         & 5x,'      Profile (average)=    ',es11.4, /,&
         & 5x,'      --------------------- ',        /,&
         & 5x,'      Total volume=         ',es11.4,/,&
         & 5x,'      Ave. element volume=  ',es11.4,/,&
         & 5x,'      Min. element volume=  ',es11.4,' (element ',a,')',/,&
         & 5x,'      Max. element volume=  ',es11.4,' (element ',a,')',/,&
         & 5x,'      --------------------- ',   /,&
         & 5x,'      Total elements=       ',i11,/,&
         & 5x,'      Extension elements=   ',i11,/,&
         & 5x,'      Hole elements=        ',i11,/,&
         & 5x,'      Elements used:        ',a, /,&
         & 5x,'      Gauss points (rule)=  ',a, /,&
         & 5x,'      --------------------- ',   /,&
         & 5x,'      Total boundaries=     ',i11,/,&
         & 5x,'      Extension boundaries= ',i11,/,&
         & 5x,'      Hole boundaries=      ',i11,/,&
         & 5x,'      Boundaries used:      ',a, /,&
         & 5x,'      Gauss points (rule)=  ',a, /,&
         & 5x,'      --------------------- ',   /,&
         & 5x,'      Total nodes=          ',i11,/,&
         & 5x,'      Hole nodes=           ',i11)

1701 format(&
         & 5x,'      --------------------- ',   /,&
         & 5x,'      Total edges=          ',i11)

1800 format(//,&
         & 5x,'   --------------------------',/&
         & '*',4x,'|- SUMMARY OF COMPUTING TIMES',/,&
         & 5x,'   --------------------------',/,&
         & /,&
         & /,&
         & 5x,'   TOTAL CPU TIME:              ',F10.2,/,&
         & 5x,'   STARTING OPERATIONS:         ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Geometry reading:          ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Set Reading:               ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Boundary cond. reading:    ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Fields reading:            ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Mesh partitioning:         ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Mesh multiplication:       ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Domain construction:       ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Groups:                  ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Halos:                   ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Element search:          ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Coupling:                ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Output mesh:             ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Additional arrays:         ',F10.2,' (',F6.2,' % )')
1900 format(/,&
         & 5x,'   ',a6,' MODULE:               ',F10.2,' (',F6.2,' % )')
1901 format(&
         & 5x,'     Parallel efficiency:       ',F10.2,' %')
9900 format(&    
         & 5x,'     ',a,/,&
         & 5x,'       Average:                 ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Maximum:                 ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Average load balance:    ',F10.2,/,&
         & 5x,'       Average element cost:    ',F10.2,' ',a,'s',/,&
         & 5x,'       ---                      ',/,&
         & 5x,'       Theoretical average:     ',F10.2,/,&
         & 5x,'       Theoretical maximum:     ',F10.2,/,&
         & 5x,'       Theoretical load bal:    ',F10.2,/,&
         & 5x,'       Maximum variability:     ',F10.2,' % ')
10000 format(&   
         & 5x,'     ',a,/,&
         & 5x,'       Average:                 ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Maximum:                 ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Load balance:            ',F10.2,/,&
         & 5x,'       Average entity cost:     ',F10.2,' ',a,'s')
10001 format(&   
         & 5x,'     ',a,/,&
         & 5x,'       Average:                 ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Maximum:                 ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'       Load balance:            ',F10.2)
10002 format(&   
         & 5x,'     ',a,F10.2,' (',F6.2,' % )')
2000 format(&
         & 5x,'   OUTPUT OPERATIONS:           ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   OTHER OPERATIONS:            ',F10.2,' (',F6.2,' % )')
2001 format(&
         & 5x,'   OTHER OPERATIONS:            ',F10.2,' (',F6.2,' % )')
2100 format(//,&
         & 5x,'   --------------------------',/&
         & '*',4x,'|- SUMMARY OF REQUIRED MEMORY',/,&
         & 5x,'   --------------------------',/,&
         & /)
2102 format(&
         & 5x,'   |- MASTER COUNT:',/)
2103 format(/,&
         & 5x,'   |- SLAVE MAX COUNT:',/)
2101 format(&
         & 5x,'      Total memory:                ',f7.2,1x,a6,/,&
         & 5x,'      Construction of the domain:  ',f7.2,1x,a6,/,&
         & 5x,'      Element search strategy:     ',f7.2,1x,a6 )
2200 format(&
         & 5x,'      ',a6,' module:               ',f7.2,1x,a6 )

2400 format(&
         & 5x,'      Solvers:                     ',f7.2,1x,a6 )
2500 format(/,&
         & 5x, '   ',a,/,&
         & 5x, '|- ',a,/,&
         & 5x, '   ',a,/)
2800 format(/,&
         & 5x,'|- ',a6,' IS STATIONARY AT TIME STEP ',i5)
    !
    ! ----------------------------------------------- COMPUTING TIMES FORMATS
    !
2900 format(//,&
         & 5x,'   --------------------------',/&
         & 5x,'|- SUMMARY OF COMPUTING TIMES',/,&
         & 5x,'   --------------------------',/,&
         & //,&
         & 5x,'   TOTAL CPU TIME:            ',es10.3)
3000 format(&
         & 5x,'     ',a25,es10.3,' (',f6.2,' % )')
3100 format(&
         & 5x,'     ',a25,es10.3,' (',f6.2,' % )')
    !
    ! ---------------------------------------------- COMPUTING TIMES FORMATS
    !
3200 format(&
         & '# Iterations  = ',3(1x,i9),/,&
         & '# Time        = ',e16.8E3)
3300 format(&
         & '# NUMVARIABLES : ',i7,/,&
         & '# NUMSETS      : ',i7,/,&
         & '# START')
3400 format(&
         & '# ALYA ',a,' set results',/,&
         & '#',/,&
         & '# HEADER')
3500 format('# ',a,' , Column : ',i3)

3600 format(//,&
         & 5x,'   -----------------',/&
         & 5x,'|- SUMMARY OF SOLVER',/,&
         & 5x,'   -----------------',/,&
         & /)
3700 format(&
         & 5x,'   ',a,':',/,&
         & 5x,'   ',a,/,&
         & /,&
         & 5x,'   SOLVER:               ',a,/,&
         & 5x,'   # SOLVES=             ',i8)
3701 format(&
         & 5x,'   MIN. # ITERATIONS=    ',i8,/,&
         & 5x,'   MAX. # ITERATIONS=    ',i8,/,&
         & 5x,'   AVE. # ITERATIONS=    ',i8,/,&
         & 5x,'   TOT. # ITERATIONS=    ',i8,/,&
         & 5x,'   TOT. CPU TIME=        ',es10.3)
3702 format(&
         & 5x,'   ',a22,i8)
3703 format(&
         & 5x,'   ',a22,a)
3704 format(&
         & 5x,'   ',a22,e12.6)
3705 format(//,&
         & 5x,'   SUMMARY OF SOLVER TIMES:',/,&
         & /,&
         & 5x,'   TOTAL SOLVER CPU TIME (AVERAGE): ',es10.3,/,&
         & 5x,'     Initial operations:            ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'     Preconditioning:               ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'     Coarse system solver (DCG):    ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'     Orthogonalization (GMRES):     ',es10.3,' (',F6.2,' % )')
3708 format(&
         & 5x,'     Matrix copy and exchange       ',es10.3,' (',F6.2,' % ) (Full row format)')
3709 format(&
         & 5x,'     Others:                        ',es10.3,' (',F6.2,' % )')

3707 format(//,&
         & 5x,'   |- SpMV:                    ',i9,       /,/,&
         & 5x,'      Computation:             '           /,&
         & 5x,'         Average CPU time:     ',es10.3,' ',/,&
         & 5x,'         Maximum CPU time:     ',es10.3,' ',/,&
         & 5x,'         Load balance:         ',F10.2     /,&
         & 5x,'      Point-to-point comm.:    '           /,&
         & 5x,'         Average CPU time:     ',es10.3,' ',/,&
         & 5x,'         Maximum CPU time:     ',es10.3,' ',/,&
         & 5x,'         Load balance:         ',F10.2     /,&
         & 5x,'      Overall:                 '           /,&
         & 5x,'         Average CPU time:     ',es10.3,' ',/,&
         & 5x,'         Maximum CPU time:     ',es10.3,' ',/,&
         & 5x,'         Load balance:         ',F10.2     ,/,&
         & 5x,'                               '           ,/,&
         & 5x,'      Flops:                   '           ,/,&
         & 5x,'         Average:              ',es10.3,' ',a,'flop/s',/,&
         & 5x,'         Maximum:              ',es10.3,' ',a,'flop/s')

3717 format(//,&
         & 5x,'   |- Dot product:             ',i9,       /,/,&
         & 5x,'      Computation:             '           /,&
         & 5x,'         Average CPU time:     ',es10.3,' ',/,&
         & 5x,'         Maximum CPU time:     ',es10.3,' ',/,&
         & 5x,'         Load balance:         ',F10.2     /,&
         & 5x,'      Point-to-point comm.:    '           /,&
         & 5x,'         Average CPU time:     ',es10.3,' ',/,&
         & 5x,'         Maximum CPU time:     ',es10.3,' ',/,&
         & 5x,'         Load balance:         ',F10.2     /,&
         & 5x,'      Overall:                 '           /,&
         & 5x,'         Average CPU time:     ',es10.3,' ',/,&
         & 5x,'         Maximum CPU time:     ',es10.3,' ',/,&
         & 5x,'         Load balance:         ',F10.2     ,/,&
         & 5x,'                               '           ,/,&
         & 5x,'      Flops:                   '           ,/,&
         & 5x,'         Average:              ',es10.3,' ',a,'flop/s',/,&
         & 5x,'         Maximum:              ',es10.3,' ',a,'flop/s')

3706 format(//,&
         & 5x,'   |- TOTAL SOLVER CPU TIME:      ',es10.3,/,/,&
         & 5x,'      Invert Aii:                 ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'      Compute preconditioner:     ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'      Factorize preconditioner:   ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'      Solver initialization:      ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'      Solver main loop. q = A p:  ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'      Solver main loop. descent:  ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'      Solver main loop. L z = q:  ',es10.3,' (',F6.2,' % )',/,&
         & 5x,'      Solver main loop. others:   ',es10.3,' (',F6.2,' % )')

3800 format(//,&
         & 5x,'   --------------------',/&
         & '*',4x,'|- SUMMARY OF PARTITION',/,&
         & 5x,'   --------------------',//,&
         & 5x,'   # SUBDOMAINS=         ',i11,/,&
         & 5x,'   TOTAL # NODES=        ',i11,' (INTERIOR+REPEATED BOUNDARY)',//,&
         & 5x,'   MIN. # NPOIN=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   MAX. # NPOIN=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   AVE. # NPOIN=         ',i11,/,&
         & 5x,'   MIN. # NELEM=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   MAX. # NELEM=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   AVE. # NELEM=         ',i11,/,&
         & 5x,'   LB     NELEM=         ',F11.3,/,&
         & 5x,'   MIN. # NELEM(WEIGHT)= ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   MAX. # NELEM(WEIGHT)= ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   AVE. # NELEM(WEIGHT)= ',i11,/,&
         & 5x,'   LB     NELEM(WEIGHT)= ',F11.3,/,&
         & 5x,'   MIN. # NBOUN=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   MAX. # NBOUN=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   AVE. # NBOUN=         ',i11,/,&
         & 5x,'   MIN. # NZDOM=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   MAX. # NZDOM=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   AVE. # NZDOM=         ',i11,//,&
         & 5x,'   MIN. # NNEIG=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   MAX. # NNEIG=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   AVE. # NNEIG=         ',i11,/,&
         & 5x,'   MIN. # NBBOU=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   MAX. # NBBOU=         ',i11,' (SUBDOMAIN ',i5,')',/,&
         & 5x,'   AVE. # NBBOU=         ',i11)

3903 format(&
         & '#',/,&
         & '# -----------------',/,&
         & '# ALGEBRAIC SOLVER:',/,&
         & '# -----------------',/,&
         & '# ',/,&
         & '# TYPE OF SOLVER        : ',a,/,&
         & '# --------------')
3904 format(&
         & '# RELAXATION            = ',es10.3)
3905 format(&
         & '# KRYLOV DIMENSION      = ',i10)
3950 format(&
         & '# ORTHOGONALIZATION     = ',a)
3906 format(&
         & '# MINIMUM TOLERANCE     = ',e10.4)
3918 format(&
         & '# ADAPTIVE RESIDUAL     : OFF')
3919 format(&
         & '# ADAPTIVE RESIDUAL     : ON, RATIO=',e10.4)
3907 format(&
         & '# MAX. ITERATIONS       = ',i10)
3908 format(&
         & '# NUMBER OF GROUPS      = ',i10)
3909 format(&
         & '# TOLERANCE OF LINELETS = ',es10.3,/,&
         & '# NUMBER OF LINELETS    = ',i10,/,&
         & '# CROSSED BY BOUNDARIES = ',i10,/,&
         & '# LINELET NODES         = ',i10,' %')
3910 format(&
         & '# NON-ZERO COEFF.       = ',i10,' (R*',i1,')',', MEMORY= ',f7.2,' ',a)
3915 format(&
         & '# MIN NON-ZERO COEFF.   = ',i10,' (R*',i1,')',', MEMORY= ',f7.2,' ',a)
3916 format(&
         & '# MAX NON-ZERO COEFF.   = ',i10,' (R*',i1,')',', MEMORY= ',f7.2,' ',a)
3917 format(&
         & '# AVE NON-ZERO COEFF.   = ',i10,' (R*',i1,')',', MEMORY= ',f7.2,' ',a)
3935 format(&
         & '# FORMAT                : ',a)
3911 format(&
         & '# COARSE SOLVER         : ',a)
3912 format(&
         & '# COARSE RHS COMMUNIC.  : ',a)
3913 format(&
         & '#',/,&
         & '# COARSE CORRECTION     : ',/,&
         & '# -----------------')
3914 format(&
         & '# COARSE SPACE SIZE     = ',i10)
3929 format(&
         & '#',/,&
         & '# THIS IS A TIME RESTART RUN: CONTINUING... ',/,'#')
3920 format(&
         & '#',/,&
         & '# PRECONDITIONER        : ',a,/,&
         & '# --------------')
3921 format(&
         & '# COARSE SYSTEM SOLVER  : ',a)
3922 format(&
         & '# NUMBER OF GROUPS      = ',i10)
3923 format(&
         & '# PRECONDITIONER        = ',a)
3928 format(&
         & '# BLOCK GAUSS-SEIDEL    : ',a)
3924 format(&
         & '# SMOOTHER              = ',a)
3925 format(&
         & '# SMOOTHER ITERATIONS   = ',i10)
3926 format(&
         & '# AII SOLVER            = ',a)
3927 format(&
         & '# ABB SOLVER            = ',a)
3931 format(&
         & '# LEFT                  : ',a)
3933 format(&
         & '# RIGHT                 : ',a)
3934 format(&
         & '# LEFT/RIGHT            : ',a)
3932 format(&
         & '# DIRECT SOLVER         : ',a)
3930 format(&
         & '#',/,&
         & '# --|  Columns displayed:' ,/,&
         & '# --|  1. Number of iterations   2. Initial prec. residual   ',/,&
         & '# --|  3. Final prec. residual   4. Initial residual         ',/,&
         & '# --|  5. Final residual         6. Convergence rate         ',/,&
         & '# --|  7. CPU time               8. Orthogonality            ',/,&
         & '# --|  9. RHS norm ||b||        10. Precon. upd (1=yes,0=no) ',/,& 
         & '# --| 11. Cond. number kappa    12. Min. eigenvalue          ',/,& 
         & '# --| 13. Max. eigenvalue                                    ',/,&
         & '# --|                                                        ',/,&
         & '# --| -999 Means was not calculated                          ',/,&
         & '# --| -888 Means was not converged                           ',/,&
         & '# --|                                                        ',/,&
         & '#                                                            ',/,&
         & '#   ','  1','                 2','                 3','                 4', &
         &              '                 5','                 6','                 7', &
         &              '                 8','                 9','            10', &
         &              '                11','                12','                13',/,&
         & '#')

3940 format(&
         & '#',/,&
         & '# --| Columns displayed:' ,/,&
         & '# --| 1. Number of iter.        2. Preconditioned residual ',/,&
         & '# --| 3. Residual               4. Time of current iter.   ',/,&
         & '# --| 5. Time from first iter.  6. Orthogonality           ',/,&
         & '# ','    1','                2','                3','                4',&
         &              '                5','                6',/,&
         & '#')
3941 format(&
         & '#',/,&
         & '# PARALLELIZATION       : ',/,&
         & '# ---------------')
3942 format(&
         & '# SpMV MPI MATRIX       : ',a)
3936 format(&
         & '# SpMV MPI COMMUNICAT.  : ',a)
3943 format(&
         & '# SpMV OMP              : ',a,'IN ',i7, 'SUBDOMAINS')
3944 format(&
         & '#    SpMV OMP CHUNK SIZE : ',i7,' (MIN)')
3945 format(&
         & '#    SpMV OMP CHUNK SIZE : ',i7,' (MAX)')
3946 format(&
         & '#    SpMV OMP CHUNK SIZE : ',i7,' (AVE)')
4000 format(//,&
         & 5x,'   -------------------',/&
         & 5x,'|- SUMMARY OF REVISION',/,&
         & 5x,'   -------------------',/,&
         & /,&
         & /)
4001 format(5x,'   ',a)
4100 format(//,&
         & 5x,'   --------------------------',/&
         & 5x,'|- SUMMARY OF RUN ENVIRONMENT',/,&
         & 5x,'   --------------------------',/,&
         & /)
4101 format(8x,'   ',a)
4102 format(/,&
         & 5x,'   |- USER AND SYSTEM:',/)
4103 format(/,&
         & 5x,'   |- COMPILERS:',/)

4200 format(//,&
         & 5x,'   ------------------------',/&
         & '*',4x,'|- BOUNDARY CONDITION CODES',/,&
         & 5x,'   ------------------------',/,&
         & /)
4301 format(&
         & 5x,'   CODE=            ',(i4))
4302 format(&
         & 5x,'   CODE=            ',i4,'  & ',i4)
4303 format(&
         & 5x,'   CODE=            ',i4,'  & ',i4,'  & ',i4)
    !4301 format(&
    !       & 5x,'   NUMBER OF NODES= ',i10)

4400 format(//,&
         & 5x,'   ---------------------------',/  &
         & 5x,'|- GEOMETRICAL CONDITION NODES',/, &
         & 5x,'   ---------------------------',//,&
         & 5x,'   PRESCRIBED=       ',i10,/,&
         & 5x,'   FREESTREAM=       ',i10,/,&
         & 5x,'   WALL LAW=         ',i10,/,&
         & 5x,'   WALL LAW (LINE)=  ',i10,/,&
         & 5x,'   SYMMETRY=         ',i10,/,&
         & 5x,'   SYMMETRY (LINE)=  ',i10,/,&
         & 5x,'   NO SLIP=          ',i10,/,&
         & 5x,'   FREE/OUTFLOW=     ',i10)

4500 format(//,&
         & 5x,'   BAD EDGES:        ',i10,/,&
         & 5x,'   LAST EDGE FOUND:  ',i10,'-',i10,/,&
         & 5x,'   LAST BOUNDARY=    ',i10)
4600 format(&
         & '# Iteration  = ',1(1x,i6),/,&
         & '# Time       = ',e16.8E3)
4900 format(//,&
         & 5x,'   -------------------',/  &
         & '*',4x,'|- MESH MULTIPLICATION',/, &
         & 5x,'   -------------------',//,&
         & 5x,'   EDGE TABLE:               ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   FACE TABLE:               ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   ELEMENT DIVISION:         ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   PARALL, OWN NODES:        ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   PARALL, COMMON EDGES:     ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   PARALL, COMMON FACES:     ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   PARALL, BOUNDARIES:       ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   PARALL, LOC. RENUMBERING: ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'   PARALL, GLO. NUMBERING:   ',F10.2,' (',F6.2,' % )',/)
5000 format(//,&
         & 5x,'   --------------',/&
         & '*',4x,'|- WITNESS POINTS',/,&
         & 5x,'   --------------',/)
5100 format(5x,i9)
5101 format(//,&
         & 9x,'SOME WITNESS POINTS DO NOT HAVE HOST ELEMENTS:')
5200 format(//,&
         & 9x,'ALL WITNESS POINTS HAVE FOUND HOST ELEMENTS')
5300 format('# --| ALYA Lagrangian postprocess  ' ,/,&
         & '# --| Columns displayed:' ,/,&
         & '# --|  1. Current time          2. Particle number        3. x            ',/,&
         & '# --|  4. y                     5. z                      6. vx           ',/,&
         & '# --|  7. vy                    8. vz                     9. Type number  ',/,&
         & '# --| 10. Subdomain            11. ax                    12. ay           ',/,&
         & '# --| 13. az                   14. Cd                    15. Stokes nb.   ',/,&
         & '# --| 16. vfx                  17. vfy                   18. vfz          ',/,&
         & '# --| 19. adx                  20. ady                   21. adz          ',/,&
         & '# --| 22. aex                  23. aey                   24. aez          ',/,&
         & '# --| 25. agx                  26. agy                   27. agz          ',/,&
         & '# --| 28+ physical properties')
5301 format(&
         & '# ','         1','            2','                 3',                                 &
         &      '                 4','                 5','                 6','                 7',       &
         &      '                 8','        9','       10','                11','                12',&
         &      '                13','                14','                15','                16',       &
         &      '                17')

5302 format(&
         & '# ALYA Lagrangian postprocess  ',/,&
         & '#',/,&
         & '# HEADER')
5400 format('# --| ALYA Lagrangian convergence  ' ,/,&
         & '# --| Columns displayed:' ,/,&
         & '# --|  1. Current time          2. Accumulated particles  3. Existing particles',/,&
         & '# --|  4. Deposited particles   5. Particles with dt=0    6. CPU time max      ',/,&
         & '# --|  7. CPU time avg          8. Load balance           9. Max particles sent',/,&
         & '#---| 10. Max particles recv   11. Comm. loops           ',/,&
         & '# ','         1','          2','        3','        4','        5','             6',&
         &      '             7','             8','        9','       10','       11')



5500 format('# --| ALYA Lagrangian convergence with thermodynamic data  ' ,/,&
         & '# --| Columns displayed:' ,/,&
         & '# --|  1. Current time          2. Accumulated particles  3. Existing particles',/,&
         & '# --|  4. Deposited particles   5. Particles with dt=0    6. CPU time max      ',/,&
         & '# --|  7. CPU time avg          8. Load balance           9. Max particles sent',/,&
         & '#---| 10. Max particles recv   11. Comm. loops           12. Total mass of pts ',/,&
         & '#---| 13. Total enth of pts    14. Total mom of pts x    15. Total mom of pts y',/,&
         & '#---| 16. Total mom of pts z   17. Total kine of pts     18. Bounding box x min',/,&
         & '#---| 19. Bounding box y min   20. Bounding box z min    21. Bounding box x max',/,&
         & '#---| 22. Bounding box y max   23. Bounding box z max                          ',/,&
         & '# ','         1','          2','        3','        4','        5','             6',&
         &      '             7','             8','        9','       10','       11','            12',&
         &      '            13','            14','            15','            16','            17',&
         &      '            18','            19','            20','            21','            22',&
         &      '            23')


5600 format(//,&
         & 5x,'   ---------------------------',/&
         & 5x,'|- OTHER TASKS COMPUTING TIMES',/,&
         & 5x,'   ---------------------------')

5700 format('Current time,Particle number,Particle type,x,y,z,wall(-2)/outflow(-4),boundary set')
    !5700 format('# --| ALYA Lagrangian deposition  ' ,/,&
    !       & '# --| Columns displayed:' ,/,&
    !       & '# --|  1. Current time          2. Particle number        3. Particle type',/,&
    !       & '# --|  3. x                     4. y                      5. z')
    !5701 format(&
    !       & '# ','         1','        2','        3',                                 &
    !       &      '                4','             5')
5800 format(//,&
         & 5x,'   -----------------------------',/  &
         & 5x,'|- PROJECTION ON SUPPORT SURFACE',/, &
         & 5x,'   -----------------------------',//,&
         & 5x,'   MIN DISLACEMENT = ',es10.3,' % ',/,&
         & 5x,'   MAX DISLACEMENT = ',es10.3,' % ',/,&
         & 5x,'   AVE DISLACEMENT = ',es10.3,' % ')

5900 format(//,&
         & 5x,'   -----------------------------',/  &
         & 5x,'|- WAKE MODEL INFORMATION'       ,/, &
         & 5x,'   -----------------------------',//)
6000 format(&
         & /,&
         & 5x,'   WIND TURBINE:               ',i4,/,&
         & 5x,'     A_disk=                   ',1(1x,f15.3),/,&
         & 5x,'     d_ref=                    ',1(1x,f15.7),/,&
         & 5x,'     n_disk=                   ',3(1x,f15.7),/,&
         & 5x,'     x_disk=                   ',3(1x,f15.7),/,&
         & 5x,'     x_ref=                    ',3(1x,f15.7),/,&
         & 5x,'     V_disk=                   ',1(1x,f15.3),/,&
         & 5x,'     d_disk=                   ',1(1x,f15.7))

6100 format(&
         & 8x,a31,i9)
6200 format(&
         & 8x,a31,f6.2)
6300 format(//,&
         & 5x,'   -------------------',/&
         & 5x,'|- ZONE-WISE PARTITION',/,&
         & 5x,'   -------------------',/)
6400 format(&
         & 5x,'   ZONE, # SUBDOMAINS    ',i3,',',i8)
6500 format(&
         & /,&
         & 5x,'   STEP:                 ',i4,/,&
         & 5x,'   WIND TURBINE:         ',i4,/,&
         & 5x,'   Ct =                  ',1(1x,f10.3),/,&
         & 5x,'   Cp =                  ',1(1x,f10.3),/,&
         & 5x,'   Veinf =               ',1(1x,f10.3),/,&
         & 5x,'   Hub Velocity  =       ',1(1x,f10.3),/,&
         & 5x,'   Averaged Veloc =      ',1(1x,f10.3),/,&
         & 5x,'   Power(kW) =           ',1(1x,f10.0),/,&
         & 5x,'   Yaw angle(degrees) =  ',1(1x,f10.1))


6600 format(//,&
         & 5x,'   ------------------------',/,&
         & '*',4x,'|- INTER COLOR COMMUNICATOR',/,&
         & 5x,'   ------------------------',/)
6601 format(&
         & /,&
         & 5x,'   TOTAL NUMBER OF CODES=      ',i4,/,&
         & 5x,'   TOTAL NUMBER OF ZONES=      ',i4,/,&
         & 5x,'   TOTAL NUMBER OF SUBDS=      ',i4,/,&
         & /,&
         & 5x,'   NUMBER OF ZONES IN MY CODE= ',i4,/,&
         & 5x,'   NUMBER OF SUBDS IN MY CODE= ',i4,/,&
         & 5x,'                               ',/,&
         & 5x,'   LIST OF COMMUNICATORS:      ',/)
6700 format(&
         & 5x,'   ',a,'(',i3,',',i3,'), ',a,'(',i3,',',i3,'), ',a,'(',i3,',',i3,')')

6800 format(//,&
         & 5x,'   ------------------------',/,&
         & '*',4x,'|- PARTITIONS BIN STRUCTURE',/,&
         & 5x,'   ------------------------',/,&
         & /,&
         & 5x,'   BIN DIMENSIONS=             ',i5,i5,i5)
6900 format(/,&
         & 5x,'   MAXIMUM OVER SLAVES:        ',/&
         & 5x,'   --------------------        ',//&
         & 5x,'   TOTAL MEMORY:               ',f7.2,1x,a6,/&
         & 5x,'   MAXIMUM MEMORY REQUIRED:    ',f7.2,1x,a6)

7000 format(//,&
         & 5x,'   -----------------------------',/,&
         & 5x,'|- LOAD BALANCE & COMMUNICATIONS',/,&
         & 5x,'   -----------------------------',/,&
         & /,&
         & 5x,'    # % elements % boundaries % neighbors')
7100 format(&
         & 5x,3x,3(1x,es10.3))
7200 format(//,&
         & 5x,    '   ----------------------',/,&
         & '*',4x,'|- HYBRID PARALLELIZATION',/,&
         & 5x,    '   ----------------------',/)
7206 format(&
         & 5x,'   OFF')
7205 format(&
         & 5x,'   STRATEGY:                   ',/,&
         & 5x,'   LOOPS WITH RACE CONDITION:  ',a,/,&
         & 5x,'   LOOPS WITHOUT RACE:         ',a,/,&
         & 5x,'                               ',/,&
         & 5x,'   DYNAMIC LOAD BALANCE (DLB): ',a,/,&
         & 5x,'                               ',/,&
         & 5x,'   NUMBER OF THREADS=          ',i7,/,&
         & 5x,'   NUMBER OF BLOCKS=           ',i7,/,&
         & 5x,'   GRANULARITY=                ',i7,/,&
         & 5x,'   AVERAGE ELEMENT CHUNK=      ',i7,/,&
         & 5x,'   AVERAGE BOUNDARY CHUNK=     ',i7,/,&
         & 5x,'   AVERAGE NODE CHUNK=         ',i7)
7202 format(&
         & 5x,'   AVERAGE # OF COLORS=        ',i7)
7204 format(&
         & 5x,'   AVERAGE # OF SUBDOMAINS=    ',i7,/,&
         & 5x,'   AVERAGE # OF DEPENDENCIES=  ',i7)
7203 format(&
         & 5x,'    #  Color  number_of_elements')

7201 format(&
         & 5x,3x,2(1x,i9))
7210 format(//,&
         & 5x,    '   ------------',/,&
         & '*',4x,'|- CO-EXECUTION',/,&
         & 5x,    '   ------------',/)
7211 format(&
         & 5x,'   # OF CPU PROCESSES=         ',i7,/,&
         & 5x,'   # OF GPU PROCESSES=         ',i7)

7300 format(//,&
         & 5x,    '   --------------------',/,&
         & '*',4x,'|- OPTIMIZATION OPTIONS',/,&
         & 5x,    '   --------------------',/,&
         & //,&
         & 5x,'   MPI3:                       ',a3,/&
         & 5x,'   BLAS:                       ',a3)
7305 format(&
         & 5x,'   VECTOR SIZE=                ',i8,' (VARIABLE)')
7307 format(&
         & 5x,'   VECTOR SIZE=                ',i3,' (PARAMETER)')
7306 format(&
         & 5x,'   VECTOR SIZE:                ','ADJUSTED TO CHUNK')
7301 format(&
         & 5x,'   LOC. TO GLOB. COEF IN CSR:  OFF')
7302 format(&
         & 5x,'   LOC. TO GLOB. COEF IN CSR:  ON, AVERAGE MEMORY= ',f7.2,' ',a)
7303 format(&
         & 5x,'   ELEMENT DATA BASE:          OFF')
7304 format(&
         & 5x,'   ELEMENT DATA BASE:          ON, AVERAGE MEMORY= ',f7.2,' ',a)
7308 format(&
         & 5x,'   GPU/CPU MPI PROCESSES=      ',i6,', ',i6)
7309 format(&
         & 5x,'   NODE RENUMBERING:           ',a)
7310 format(&
         & 5x,'   ELEMENT RENUMBERING:        ',a)

7400 format(//,&
         & 5x,    '   -------------------------',/&
         & '*',4x,'|- SUMMARY OF DIRECT SOLVERS ',/,&
         & 5x,    '   -------------------------')

7500 format(//,&
         & 5x,'    |- MODULE:                     ',a,/,&
         & 5x,'       EQUATION:                   ',a,/,&
         & 5x,'       PROBLEM:                    ',a,/,&
         & 5x,'       DIRECT SOLVER:              ',a,/,&
         & 5x,'       NUMBER THREADS:             ',i6,/,&
         & 5x,'       INITIALIZATION STEP:        ',/,&
         & 5X,'          Number:                  ',i9,/,&
         & 5x,'          Average CPU time:        ',F10.2,/,&
         & 5x,'          Maximum CPU time:        ',F10.2,/,&
         & 5x,'          Load balance:            ',F10.2,/,&
         & 5x,'       FACTORIZATION STEP:         ',/,&
         & 5X,'          Number:                  ',i9,/,&
         & 5x,'          Average CPU time:        ',F10.2,/,&
         & 5x,'          Maximum CPU time:        ',F10.2,/,&
         & 5x,'          Load balance:            ',F10.2,/,&
         & 5x,'          Average fill-in:         ',F10.2,/,&
         & 5x,'          Maximum fill-in:         ',F10.2,/,&
         & 5x,'       SOLUTION STEP:              ',/,&
         & 5X,'          Number:                  ',i9,/,&
         & 5x,'          Average CPU time:        ',F10.2,/,&
         & 5x,'          Maximum CPU time:        ',F10.2,/,&
         & 5x,'          Load balance:            ',F10.2,/,&
         & 5x,'       MAX MEMORY ALLOCATED:       ',/,&
         & 5x,'          Average=                 ',f7.2,' ',a,/,&
         & 5x,'          Maximum=                 ',f7.2,' ',a)

7600 format(//,&
         & 5x,    '   -----------------------------',/&
         & '*',4x,'|- TIMES OF ALGEBRAIC OPERATIONS: AVERAGE OVER ',i5,' RUNS',/,&
         & 5x,    '   -----------------------------',/,&
         & /,&
         & /,&
         & 5x,'   |- SpMV IN CSR FORMAT (without OpenMP):   ',/,/,&
         & 5x,'      Computation:             '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Point-to-point comm.:    '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Overall:                 '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5         )
7601 format(/,&
         & 5x,'   |- SpMV IN CSR FORMAT (with OpenMP & static scheduling):',/,/,&
         & 5x,'      Computation:             '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Point-to-point comm.:    '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Overall:                 '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5         )
7604 format(/,&
         & 5x,'   |- SpMV IN CSR FORMAT (on GPU):'            ,/,/,&
         & 5x,'      Data transfer CPU->GPU:  '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Computation:             '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Overall:                 '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5         )
7602 format(/,&
         & 5x,'   |- DOT PRODUCT (without OpenMP):'        ,/,/,&
         & 5x,'      Computation:             '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Global comm.:            '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Overall:                 '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5         )
7603 format(/,&
         & 5x,'   |- DOT PRODUCT (with OpenMP):   '        ,/,/,&
         & 5x,'      Computation:             '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Global comm.:            '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Overall:                 '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5         )
7605 format(/,&
         & 5x,'   |- DOT PRODUCT (on GPU):    '          ,/,/,&
         & 5x,'      Data transfer CPU->GPU:  '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Global comm.:            '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5      ,/,&
         & 5x,'      Overall:                 '            ,/,&
         & 5x,'         Average CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:     ',F10.5,' ',a,/,&
         & 5x,'         Load balance:         ',F10.5         )

7700 format(//,&
         & 5x,    '   -----------------',/,&
         & '*',4x,'|- MESH PARTITIONING',/,&
         & 5x,    '   -----------------',/,/,&
         & 5x,'   METHOD:                     ',a)
7708 format(&
         & 5x,'   EXECUTION MODE              ',a)
7701 format(&
         & 5x,'   NUMBER OF BOXES=            ',3(1x,i7))
7702 format(&
         & 5x,'   NUMBER OF BOXES (FINE)=     ',3(1x,i7))
7703 format(&
         & 5x,'   NUMBER OF BOXES (COARSE)=   ',3(1x,i7))
7704 format(&
         & 5x,'   DIRECTION OF PARTITIONING=  ','(',e13.6,',',e13.6,',',e13.6,')')
7707 format(&
         & 5x,'   DIRECTION OF PARTITIONING=  ','(',e13.6,',',e13.6,')')
7705 format(&
         & 5x,'   GRAPH NEIGHBOR CRITERION:   ',a)
7706 format(&
         & 5x,'   WEIGHT CRITERION FOR LB:    ',a)
7709 format(&
         & 5x,'   TOTAL TIME:                ',es10.3,/,&
         & 5x,'     Redistribution from I/O= ',es10.3,' (',F6.2,' % )',/,& ! 1
         & 5x,'     Weights calculation=     ',es10.3,' (',F6.2,' % )',/,& ! 2
         & 5x,'     Partitioning=            ',es10.3,' (',F6.2,' % )',/,& ! 3
         & 5x,'     Mesh redistribution=     ',es10.3,' (',F6.2,' % )',/,& ! 4
         & 5x,'     Communication arrays=    ',es10.3,' (',F6.2,' % )')    ! 5

7800 format(//,&
         & 5x,    '   ---------------------',/,&
         & '*',4x,'|- INTEGRATION RULE TEST',/,&
         & 5x,    '   ---------------------',/,/)
7900 format(&
         & 5x,'   ELEMENT, ORDER, GAUSS, VOLUME, ERROR=    ',a,', ',i2,', ',i3,', ',e13.6,', ',e13.6,' -> ',a)

8000 format(//,&
         & 5x,    '   -----------------',/,&
         & '*',4x,'|- INTEGRATION RULES',/,&
         & 5x,    '   -----------------',/,/)
8100 format(&
         & 5x,'   ELEMENT, RULE, GAUSS POINTS, MAX. ORDER = ',a,', ',a,', ',i3,', ',i2)
8101 format(&
         & 5x,'   ELEMENT =                               = WRONG INTEGRATION RULE! CHECK INPUT DATA')
8102 format(&
         & 5x,'   ELEMENT TYPE:                 ',a,/,&
         & 5x,'      Rule:                      ',a,/,&
         & 5x,'      Number Gauss points=       ',i5,/,&
         & 5x,'      Max. polynomial order=     ',i5)
8103 format(&
         & 5x,'   ELEMENT TYPE:                 ',a,/,&
         & 5x,'       WRONG INTEGRATION RULE! CHECK INPUT DATA')

8200 format(//,&
         & 5x,    '   -----',/&
         & '*',4x,'|- HALOS',/,&
         & 5x,    '   -----',/,&
         & /,&
         & /,&
         & 5x,'   |- COMPUTATIONAL HALOS (MAX VALUES):'            ,/,/,&
         & 5x,'      Halo nodes/total nodes=       ',F10.5,' %',/,&
         & 5x,'      |neigh send - neigh recv|=    ',F10.5,' % '/,&
         & 5x,'      |size  send - size  recv|=    ',F10.5,' % '/,&
         & 5x,'      Ratio 2nd neighbors/neighbors=',F10.5,' % '/,&
         & 5x,/,&
         & 5x,'   |- GEOMETRICAL HALOS (MAX VALUES): '           ,/,/,&
         & 5x,'      Halo elements/own elements=   ',F10.5,' %',/,&
         & 5x,'      Halo nodes/own nodes=         ',F10.5,' % ')
8300 format(//,&
         & 5x,'   ---------------------------------------------------',/&
         & 5x,'|- SUMMARY OF ASSEMBLY COMPUTING TIMES (AVERAGE TIMES)',/,&
         & 5x,'   ---------------------------------------------------',/,&
         & //,&
         & 5x,'   NUMBER OF ASSEMBLIES:      ',i10,/,&
         & 5x,'   TOTAL TIME:                ',es10.3,/,&
         & 5x,'     Gather=                  ',es10.3,' (',F6.2,' % )',/,& ! 1
         & 5x,'     Element Jacobian=        ',es10.3,' (',F6.2,' % )',/,& ! 2
         & 5x,'     Properties=              ',es10.3,' (',F6.2,' % )',/,& ! 3
         & 5x,'     Residual and GP values=  ',es10.3,' (',F6.2,' % )',/,& ! 4
         & 5x,'     Subgrid scale=           ',es10.3,' (',F6.2,' % )',/,& ! 5
         & 5x,'     Stabilization parameters=',es10.3,' (',F6.2,' % )',/,& ! 6
         & 5x,'     Element matrix=          ',es10.3,' (',F6.2,' % )',/,& ! 7
         & 5x,'     Scatter in global matrix=',es10.3,' (',F6.2,' % )',/,& ! 8
         & 5x,'     Other (not in elmope)=   ',es10.3,' (',F6.2,' % )')    ! 9

8400 format(//,&
         & 5x,    '   --------------------',/&
         & '*',4x,'|- NODE CODES USED FOR ',a,/,&
         & 5x,    '   --------------------',/,&
         & /)

8500 format(//,&
         & 5x,    '   ------------------------',/&
         & '*',4x,'|- BOUNDARY CODES USED FOR ',a,/,&
         & 5x,    '   ------------------------',/,&
         & /)

8600 format(//,&
         & 5x,    '   --------',/&
         & '*',4x,'|- LOCALITY',/,&
         & 5x,    '   --------',/,&
         & /,&
         & /,&
         & 5x,'   |- Spatial of elements in assembly (number of nodes) ',/,/,&
         & 5x,'      Average:                 ',i6,/,&
         & 5x,'      Maximum:                 ',i6,/,/,&
         & 5x,'   |- Temporal of nodes in assembly (number of elements)',/,/,&
         & 5x,'      Average:                 ',i6,/,&
         & 5x,'      Maximum:                 ',i6,/,/,&
         & 5x,'   |- Spatial of nodes in SpMV (number of nodes)',/,/,&
         & 5x,'      Average:                 ',i6,/,&
         & 5x,'      Maximum:                 ',i6,/,/,&
         & 5x,'   |- Temporal of edges in assembly (number of elements)',/,/,&
         & 5x,'      Average:                 ',i6,/,&
         & 5x,'      Maximum:                 ',i6)

8700 format(//,&
         & 5x,    '   ---------------------------------',/&
         & '*',4x,'|- AUTO TUNING OF OPENMP FOR SOLVERS',/,&
         & 5x,    '   ---------------------------------',/)
8800 format(/,&
         & 5x,'   |- PROBLEM: ',a,/)
8902 format(&
         & 5x,'      Guided scheduling=  ',i7,' subdomains')
8900 format(&
         & 5x,'      Static scheduling=  ',i7,' subdomains')
8904 format(&
         & 5x,'      OpenMP desactivated=',i7,' subdomains')
8901 format(&
         & 5x,'      Dynamic scheduling= ',i7,' subdomains',/,&
         & 5x,'         Average optimum chunk size: ',i10,/,&
         & 5x,'         Minimum optimum chunk size: ',i10,/,&
         & 5x,'         Maximum optimum chunk size: ',i10)
8903 format(/,&
         & 5x,'      Timings:',/,&
         & 5x,'         Average CPU time:           ',F10.5,' ',a,/,&
         & 5x,'         Minimum CPU time:           ',F10.5,' ',a,/,&
         & 5x,'         Maximum CPU time:           ',F10.5,' ',a,/,&
         & 5x,'         Load balance:               ',F10.5)

9000 format(//,&
         & 5x,'   ---------------------------',/&
         & 5x,'|- SUMMARY OF COMPUTING TIMES:',a,/,&
         & 5x,'   ---------------------------',/,&
         & //,&
         & 5x,'   TOTAL CPU TIME:            ',es10.3)
9100 format(//,&
         & 5x,'   ----------------------',/&
         & '*',4x,'|- SUMMARY OF CONVERGENCE',/,&
         & 5x,'   ----------------------',/)

9200 format(//,&
         & 5x,    '   --------------------------------',/&
         & '*',4x,'|- AFFINITY: NUMBER OF HOSTS=',i6/,&
         & 5x,    '   --------------------------------')
9300 format(/,&
         & 5x,'   Hostname: ',a,/)
9400 format('        ')
9401 format(i6,1x)
9500 format(&
         & 5x,'   SFC PARTITIONING STATISTICS',/,&
         & 5x,'     Total number of bins=    ',es10.3,/,&
         & 5x,'     Ratio touched/total bin= ',es10.3,/,&
         & 5x,'     Max weight=              ',es10.3,/,&
         & 5x,'     Ave weight=              ',es10.3,/,&
         & 5x,'     Granularity (max/sum)=   ',es10.3)

9600 format(//,&
         & 5x    ,'   -------------------------------',/&
         & '*',4x,'|- ELEMENT WEIGHTS DURING ASSEMBLY',/,&
         & 5x,    '   -------------------------------',/,&
         & /)
9700 format(&
         & 5x,'   ',a6,' MODULE: ERROR= ',F6.2,' % ')
9701 format(&
         & 5x,'    ',a5,'= ',F6.2)

9800 format(&
         & '#',/,&
         & '# --|  Columns displayed:' ,/,&
         & '# --|  1. Time step                 2. Global iteration        3. Current time          4. Overall time',/,&
         & '# --|  5. Count element assembly    6. Max element assembly    7. Ave element assembly  ',/,&
         & '# --|  8. Count boundary assembly   9. Max boundary assembly  10. Ave boundary assembly ',/,&
         & '# --| 11. Count node assembly      12. Max node assembly      13. Ave node assembly     ',/,&
         & '# --| 14. Count particle assembly  15. Max particle assembly  16. Ave particle assembly ',/,&
         & '# --| 17. Solvers',/,&
         & '# --------- --------- ------------- ------------- ------------- ------------- ------------- ',&
         &   '------------- ------------- ------------- ',&
         &   '------------- ------------- ------------- ',&
         &   '------------- ------------- ------------- ',&
         &   '------------- ',/,&               
         & '# Time step Global it Current time  Overall time  # elem. ass.  Max elem. ass Ave elem. ass ',&
         & '# boun. ass.  Max boun. ass Ave boun. ass ',&
         & '# node  ass.  Max node  ass Ave node  ass ',&
         & '# part. ass.  Max part. ass Ave part  ass ',' Solvers',/,&
         & '# --------- --------- ------------- ------------- ------------- ------------- ------------- ',&
         &   '------------- ------------- ------------- ',&
         &   '------------- ------------- ------------- ',&
         &   '------------- ------------- ------------- ',&
         &   '------------- ')
10100 format('# --| ALYA Deposition surface ' ,i6/,&
         & '# --| Columns displayed:' ,/,&
         & '# --|  1. Current time          2. Total surface          3... All types',/,&
         & '# ','             1','                2','                3')

10200 format(//,&
         & 5x,'   ------------------',/&
         & 5x,'|- SUMMARY OF MODULES',/,&
         & 5x,'   ------------------',/,&
         & /)
10201 format(8x,'   ',a)

10300 format(//,&
         & 5x,    '   -------------------------',/&
         & '*',4x,'|- SUMMARY OF REPARTITIONING',/,&
         & 5x,    '   -------------------------',/,&
         & /,&
         & /,&
         & 5x,'   NUMBER OF REPARTITIONS :     ',i5,/,&
         & 5x,'   TOTAL CPU TIME:              ',F10.2,/,&
         & 5x,'     Destructor:                ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Partitioning:              ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Domain:                    ',F10.2,' (',F6.2,' % )',/,&
         & 5x,'     Redistribution:            ',F10.2,' (',F6.2,' % )')

10400 format(//,&
         & 5x,    '   --------------------',/&
         & '*',4x,'|- PARALLEL PERFORMANCE',/,&
         & 5x,    '   --------------------',/,&
         & /,&
         & /,&
         & 5x,    '   Definitions taken from PoP project: https://pop-coe.eu/node/69 ',/,&
         & 5x,    '                                                                  ',/,&
         & 5x,    '       t_i^w      t_i^MPI                                         ',/,&
         & 5x,    '   <------------><------->                                        ',/,&
         & 5x,    '                                                                  ',/,&
         & 5x,    '   +---------------------+  -+                                    ',/,&
         & 5x,    '   |\\\\\\\\\\\\\        +   |                                    ',/,&
         & 5x,    '   +---------------------+   |                                    ',/,&
         & 5x,    '   |\\\\\\\\\\\\\\\\\    +   | P processors                       ',/,&
         & 5x,    '   +---------------------+   |                                    ',/,&
         & 5x,    '   |\\\\\\\              +   |                                    ',/,&
         & 5x,    '   +---------------------+  -+                                    ',/,&
         & 5x,    '                                                                  ',/,&
         & 5x,    '   <---------------------> te                                     ',/,&
         & 5x,    '                                                                  ',/,&
         & 5x,    '   t_i^w ..... working time of process i                          ',/,&
         & 5x,    '   t_i^MPI ... MPI communication time of process i                ',/,&
         & 5x,    '   t_e ....... Elpased time (same for all)                        ',/,&
         & 5x,    '                                                                  ',/,&
         & 5x,    '   t_ave^w = sum_i t_i^w / P (average working time)               ',/,&
         & 5x,    '   t_max^w = max_i t_i^w     (max working time)                   ',/,&
         & 5x,    '                                                                  ',/,&
         & 5x,    '                                                                  ',/,&
         & 5x,    '   Load balance:             LB = t_ave^w /t_max^w                ',/,&
         & 5x,    '   Communication efficiency: CE = t_max^w /t_e                    ',/,&
         & 5x,    '   Parallel efficiency:      PE = LB * CE = t_ave^w / te          ',/,&
         & /,&
         & /)

10401 format(&
         & 5x,'  ',a14,/,&
         & 5x,'     LB              = ',F10.2,/,&
         & 5x,'     CE              = ',F10.2,/,&
         & 5x,'     PE              = ',F10.2,/,&
         & 5x,'     Ave. time comp. = ',F10.2,/,&
         & 5x,'     Max. time comp. = ',F10.2,/,&
         & 5x,'     Ave. time MPI   = ',F10.2,/,&
         & 5x,'     Max. time MPI   = ',F10.2)

10600 format(&
         & 5x,'   |- MATERIAL DATA:',/)
10700 format(&
         & 5x,'      Material= ',i4,', volume= ',es11.4)
10800 format(&
         & 5x,'   ','PROPERTIES UPDATE:           ',F10.2,' (',F6.2,' % )')

10900 format(/,& 
         & 5x,'     Details of Doiter:',/)

11000 format(//,&
         & 5x,'   --------------------------',/&
         & '*',4x,'|- SUMMARY OF COMPUTING TIMES',/,&
         & 5x,'   --------------------------',/,&
         & /,&
         & /,&
         & 5x,'   TOTAL CPU TIME:              ',F10.2)

  end subroutine outfor

end module mod_outfor
!> @}
