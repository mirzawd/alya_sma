!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_cvgunk(itask)

  !-----------------------------------------------------------------------
  !
  ! This routine performs several convergence checks for the temperature
  !
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_solver,           only : solver_array_residual_norm
  use mod_communications,   only : PAR_SUM
  use mod_outfor,           only : outfor
  use mod_array_operations, only : array_operations_residual_norm
  use mod_array_operations, only : array_operations_min_max
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0, jpass=0
  real(rp),    save       :: cpuit_tem=0.0_rp
  real(rp)                :: ritem,time1,ritl2

  select case ( itask )

  case ( ITASK_ENDINN )
     !
     ! Check convergence of the inner iterations:
     ! || T(n,i,j) - T(n,i,j-1)|| / ||T(n,i,j)||
     !
     if( kfl_normc_tem == 3 ) then  
        ritem = solve(1) % resin   ! Algebraic residual 
     else
        ritl2 = array_operations_residual_norm(kfl_normc_tem,1_ip,1_ip,unkno,therm,0_ip,0_ip,1_ip)       
        ritem = ritl2              ! L2 residual
     end if
     if((ritem<cotol_tem).or.(itinn(modul)>=miinn_tem)) kfl_goite_tem = 0
     !
     ! Compute min and max of the temperature
     !
     if( kfl_regim_tem == 4 ) then
        call array_operations_min_max(1_ip,npoin,0_ip,tempe,temin_tem,temax_tem)
     else
        call array_operations_min_max(1_ip,nunkn_tem,0_ip,unkno,temin_tem,temax_tem)
     end if
     !
     ! Compute SGS residual
     !
     if( kfl_sgsno_tem == 1 .or. kfl_sgsti_tem == 1 ) then
        call PAR_SUM(2_ip,resgs_tem)
        if(resgs_tem(2)>0.0_rp) resgs_tem(1) = sqrt(resgs_tem(1)/resgs_tem(2))      
     end if
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        if(ipass==0) then
           time1=time1-cpu_initi
        else
           time1=time1-cpuit_tem
        end if
        if(ipass==0.and.kfl_rstar/=2) write(momod(modul)%lun_conve,100)
        write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,ritem,resgs_tem(1),&
             temin_tem,temax_tem,time1
        call cputim(cpuit_tem)
        flush(momod(modul)%lun_conve)
        ipass=1

     end if

  case ( ITASK_ENDITE )
     !
     ! Check convergence of the outer iterations:
     ! || T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||
     !
     resid_tem = array_operations_residual_norm(kfl_normc_tem,1_ip,1_ip,therm,therm,0_ip,npoin,1_ip) 

  case( ITASK_ENDSTE )
     !
     ! Check residual of the time iterations:
     ! || T(n,*,*) - T(n-1,*,*)|| / ||T(n,*,*)||
     !
     ritem = array_operations_residual_norm(kfl_normc_tem,1_ip,1_ip,therm,therm,0_ip,2_ip*npoin,1_ip)       

     if(ritem<=sstol_tem) then
        kfl_stead_tem = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if
     !
     ! Low-Mach model
     !
     if( INOTSLAVE .and. kfl_regim_tem>=3 ) then
        if( jpass == 0 ) then
           jpass = 1
           write(lun_lmach_tem,400)
        end if
        write(lun_lmach_tem,'(2x,10(2x,e12.6))') cutim,prthe(1),prthe(1)/prthe(4),dpthe,xmass_tem
     end if
  end select

  !
  ! Formats
  !
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Temperature        6. Temperaure SGS    ' ,/,& 
       &   '# --| 7. Min. temperature  8. Max. temperature   9. Elapsed CPU Time  ' ,//,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9')
101 format(4x,i9,2x,i9,2x,i9,11(2x,e12.6))
400 format('# --| ALYA Low-Mach model variables '  ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time Step         2. Therm. pres. p0   3. p0/p0^0            ',/,&
       & '# --| 4. dp/dt             5. Total mass        ',/,&
       & '# ','             1','             2','             3','             4',&
       &      '             5') 

end subroutine tem_cvgunk

