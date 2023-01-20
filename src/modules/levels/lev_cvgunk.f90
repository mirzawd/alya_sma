!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_cvgunk
  ! NAME 
  !    lev_cvgunk
  ! DESCRIPTION
  !    This routine compute the residual
  ! USES
  !    residu
  ! USED BY
  !    lev_doiter
  !***
  !-----------------------------------------------------------------------

  use      def_parame
  use      def_master
  use      def_domain
  use      def_levels
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  real(rp),    save       :: time1,time2
  real(rp)                :: rilev

  select case(itask)

  case(1)
     !
     ! Check convergence of the inner iterations:
     ! || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
     !
     call residu(kfl_normc_lev,one,one,unkno,fleve,one,one,one,1.0_rp,rilev)     
     if((rilev<cotol_lev).or.(itinn(modul)>=miinn_lev)) kfl_goite_lev = 0
     !
     ! Compute min and max of the leve amplitude
     !
     call minmax(one,npoin,zero,fleve,lemin_lev,lemax_lev)

     if(kfl_paral<=0) then
        if(ipass==0) then
           if(kfl_rstar/=2) then
              ipass=1
              write(momod(modul)%lun_conve,100)
           end if
           time1=0.0_rp
        end if
        call cputim(time2)
        time2=time2-cpuit_lev
        write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,rilev,&
             lemin_lev,lemax_lev,time2
        call cputim(cpuit_lev)
        flush(momod(modul)%lun_conve)
     end if

  case(2)
     !
     ! Check convergence of the outer iterations:
     ! || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
     !
     call residu(kfl_normc_lev,one,one,fleve(1,1),fleve(1,2),one,one,one,1.0_rp,resid_lev)

  case(3)
     !
     ! Check residual of the time iterations:
     ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
     !
     call residu(kfl_normc_lev,one,one,fleve(1,1),fleve(1,3),one,one,one,1.0_rp,rilev)
     if(rilev<=sstol_lev) then
        kfl_stead_lev = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if
  end select
  !
  ! Formats
  !
100 format('# --| Convergence '       ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       & '# --| 4. Current time      5. Leve amplitude     6. Min. leve ampl.   ' ,/,& 
       & '# --| 7. Max. leve ampl.   8. CPU time ') 
101 format(4x,i9,2x,i9,2x,i9,11(2x,e12.6))

end subroutine lev_cvgunk


