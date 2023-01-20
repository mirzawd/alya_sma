!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_measure_time.f90
!> @author  J.C. Cajas
!> @date    08/04/2015
!> @brief   Tool for time measuring
!> @details Compute the time spent in a subroutine
!> @} 
!----------------------------------------------------------------------
  subroutine PAR_MEASURE_TIME(itask, name)
    
    use def_kintyp,           only :  ip, rp
    use def_master,           only :  ini_tim, fin_tim, init_ti, fini_ti, crate_ti, cmax_ti
    use def_master,           only :  current_code, npart
    use mod_parall,           only :  PAR_MY_CODE_RANK
    use mod_parall,           only : par_memor
    use mod_memory,           only :  memory_alloca
    use mod_communications,   only :  PAR_GATHER

    implicit none
    integer(ip),   intent(in) :: itask
    character(*),  intent(in) :: name

    real(rp)                  :: dif_tim
    real(rp), pointer         :: tim_difs(:)
    integer(ip)               :: ipart

    nullify( tim_difs )
    
    if(  .not. associated(tim_difs) .and. PAR_MY_CODE_RANK == 0_ip ) call memory_alloca( par_memor, 'tim_difs', 'par_measure_time', tim_difs, npart+1_ip )

    select case( itask )

    case( 1_ip )
       !
       ! Get the initial time, to be called before enter the section of the code you are interested in
       !
       call cputim(ini_tim)
       !call system_clock(t1,rate,cmax_ti)

!        call PAR_MIN_ALL(ini_tim)

    case( 2_ip )
       !
       ! Get the final time, to be called after the section of the code you are interested in
       ! Each process calculates the time lapse and master gathers and writes the lapses 
       ! in the file fort.227+current_code
       !
       call cputim(fin_tim)
       !call system_clock(t2,rate,cmax_ti)
       dif_tim = fin_tim - ini_tim


      !       call PAR_MAX_ALL(dif_tim)
       call PAR_GATHER(dif_tim, tim_difs,'IN MY CODE')

       if( PAR_MY_CODE_RANK == 0_ip ) then 
          !call PAR_RANK_WRITE(tim_difs,name)
          open( unit=227+current_code,POSITION='append' )
          write(227+current_code,*) 'TIME SPENT IN ',name, ' FOR CODE ', current_code
          do ipart = 1_ip, npart+1_ip
            write(227+current_code,*) 'Rank= ', ipart-1_ip, tim_difs(ipart)
          end do
          flush(227+current_code)
          close( unit=227+current_code )

       end if

    case( 3_ip )
       call system_clock(init_ti,crate_ti,cmax_ti)

    case( 4_ip )
       call system_clock(fini_ti,crate_ti,cmax_ti)
       dif_tim = (fini_ti - init_ti) / crate_ti
       

!       call PAR_MAX_ALL(dif_tim)
       call PAR_GATHER(dif_tim, tim_difs,'IN MY CODE')

       if( PAR_MY_CODE_RANK == 0_ip ) then 
          !call PAR_RANK_WRITE(tim_difs,name) 
          open( unit=227+current_code,POSITION='append' )
          write(227+current_code,*) 'TIME SPENT IN ',name, ' FOR CODE ', current_code
          do ipart = 1_ip, npart+1_ip
            write(227+current_code,*) 'Rank= ', ipart-1_ip, tim_difs(ipart)
          end do
          flush(227+current_code)
          close( unit=227+current_code )

       end if
    end select 

  end subroutine PAR_MEASURE_TIME

!  subroutine PAR_RANK_WRITE(tim_difs,name)
!   
!   use def_kintyp,    only: ip,rp       
!   use def_master,    only: current_code, npart
!   
!   implicit none
!   real(rp), pointer, intent(in) :: tim_difs(:)
!   integer(ip)                   :: ipart
!   character(*),      intent(in) :: name    
!      
!   open( unit=227+current_code,access='append' )
!   write(227+current_code,*) 'TIME SPENT IN ',name, ' FOR CODE ', current_code
!   do ipart = 1_ip, npart+1_ip
!      write(227+current_code,*) 'Rank= ', ipart-1_ip, tim_difs(ipart)
!   end do
!   flush(227+current_code)
!   close( unit=227+current_code )
! end subroutine PAR_RANK_WRITE

