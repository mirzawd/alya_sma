!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_averag()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_averag
  ! NAME 
  !    tem_averag
  ! DESCRIPTION
  !    This routine averages the temperature along time
  ! USES
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod,      only : gasco 
  use mod_ker_proper,  only : ker_proper
  use mod_memory,      only : memory_alloca, memory_deallo
  use mod_solver,      only : solver_lumped_mass_system
  use mod_arrays,      only : arrays_number
  implicit none
  integer(ip)               :: ipoin, idime, dummi
  real(rp)                  :: fact0
!  real(rp)                  :: rho(2)
  real(rp), pointer         :: prope_tmp(:)

  if (cutim > avtim_tem) then

     
     if( INOTMASTER ) then


        ! AVTEM
        if( postp(1) % npp_stepi(arrays_number('AVTEM'),0) /= 0 ) then        
          ! if(kfl_regim_tem  == 3 .and. (kfl_coupl(ID_TEMPER,ID_CHEMIC) == 0)) then
          !    do ipoin=1,npoin
          !       call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
          !       avtem_tem(ipoin)=avtem_tem(ipoin)&
          !          +dtime*tempe(ipoin,1)*rho(1)
          !    end do
          ! else
              do ipoin=1,npoin
                 avtem_tem(ipoin)=avtem_tem(ipoin)&
                    +dtime*tempe(ipoin,1)
              end do
          ! end if
        end if
 
        ! FVTEM
        if( postp(1) % npp_stepi(arrays_number('FVTEM'),0) /= 0 ) then        
           if( kfl_regim_tem >= 3 ) then

             nullify ( prope_tmp )
             call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP', 'tem_averag',prope_tmp, npoin)

             call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)
             do ipoin=1,npoin
                fvtem_tem(ipoin) = fvtem_tem(ipoin) + prope_tmp(ipoin) * tempe(ipoin,1) * dtime
             end do 

             call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP', 'tem_averag',prope_tmp)
           end if
        end if
 
        ! FVTE2
        if( postp(1) % npp_stepi(arrays_number('FVTE2'),0) /= 0 ) then        
           if( kfl_regim_tem >= 3 ) then

             nullify ( prope_tmp )
             call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP', 'tem_averag',prope_tmp, npoin)

             call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)
             do ipoin=1,npoin
                fvte2_tem(ipoin) = fvte2_tem(ipoin) + prope_tmp(ipoin) &
                                    * tempe(ipoin,1) * tempe(ipoin,1) * dtime
             end do 

             call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP', 'tem_averag',prope_tmp)
           end if
        end if
      
        ! TEM**2
        if( postp(1) % npp_stepi(arrays_number('AVTE2'),0) /= 0 ) then        
          ! if(kfl_regim_tem == 3 .and. (kfl_coupl(ID_TEMPER,ID_CHEMIC) == 0)) then
          !    do ipoin=1,npoin
          !       call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
          !       avte2_tem(ipoin)=avte2_tem(ipoin)&
          !          +tempe(ipoin,1) &
          !          *tempe(ipoin,1)*dtime*rho(1)
          !    end do
          ! else
              do ipoin=1,npoin
                 avte2_tem(ipoin)=avte2_tem(ipoin)&
                    +tempe(ipoin,1) &
                    *tempe(ipoin,1)*dtime 
              end do
          ! end if
        end if

        ! VEL*TEM
        if( postp(1) % npp_stepi(arrays_number('AVTEV'),0) /= 0 ) then        
           do ipoin=1,npoin
              fact0 = tempe(ipoin,1)
              do idime =1, ndime
                 avtev_tem(idime, ipoin)=avtev_tem(idime, ipoin)&
                   +fact0*veloc(idime, ipoin,1)*dtime
              end do
           end do
        end if
        ! AVDEN
        if( postp(1) % npp_stepi(arrays_number('AVDEN'),0) /= 0 ) then
           if (kfl_regim_tem>=3) then
              fact0 = prthe(1)/gasco
              do ipoin=1,npoin
                 avden_tem(ipoin)=avden_tem(ipoin)&
                      +dtime*fact0/tempe(ipoin,1)
              end do
           end if
        end if
        ! RHO*VEL
        if( postp(1) % npp_stepi(arrays_number('FVVEL'),0) /= 0 ) then        
           if (kfl_regim_tem>=3) then
              do ipoin=1,npoin
                 fact0=prthe(1)/gasco/tempe(ipoin,1)
                 do idime =1, ndime
                    fvvel_tem(idime, ipoin)=fvvel_tem(idime, ipoin)&
                         +fact0*veloc(idime, ipoin,1)*dtime
                 end do
              end do
           end if
        end if

        ! AVRES
        if( postp(1) % npp_stepi(arrays_number('AVRES'),0) /= 0 ) then
           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP', 'tem_averag',prope_tmp, npoin)
           prope_tmp(1:npoin) =  solve(1)%reaction(1,1:npoin)
           do ipoin=1,npoin
              avres_tem(ipoin)=avres_tem(ipoin)&
                   +dtime*prope_tmp(ipoin)
           end do
           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP', 'tem_averag',prope_tmp)
        end if
        
        ! AVHSK: average heat source term from spray
        if( postp(1) % npp_stepi(arrays_number('AVHSK'),0) /= 0  ) then

           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP','tem_averag',prope_tmp,npoin)
           if (associated(heat_sink)) then
              do ipoin=1,npoin
                 prope_tmp(ipoin) = heat_sink(ipoin)
              enddo
              call solver_lumped_mass_system(1_ip,prope_tmp,EXCHANGE=.false.)
           else
              do ipoin=1,npoin
                 prope_tmp(ipoin) = 0.0_rp 
              end do
           endif

           do ipoin=1,npoin
              avhsk_tem(ipoin) = avhsk_tem(ipoin) + prope_tmp(ipoin) * dtime
           end do
           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP','tem_averag',prope_tmp)

        end if



     end if
  end if

end subroutine tem_averag
