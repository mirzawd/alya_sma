!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_updbcs()
  !-----------------------------------------------------------------------
  !****f* levels/lev_updbcs
  ! NAME 
  !    lev_updbcs
  ! DESCRIPTION
  !    This routine deals with adaptive boundary conditions
  ! USED BY
  !    lev_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_levels
  use mod_ker_space_time_function
  implicit none
  integer(ip)   :: ipoin,ibopo,idime,ifunc
  real(rp)      :: udotn,lefun,relax

  if( INOTMASTER ) then
     do ipoin = 1,npoin
        ibopo = lpoty(ipoin)
        if( ibopo /= 0 ) then
           if( abs(kfl_fixno_lev(1,ipoin)) == 8 ) then
              udotn = 0.0_rp
              do idime = 1,ndime
                 udotn = udotn + veloc(idime,ipoin,1)*exnor(idime,1,ibopo)                    
              end do
              if( udotn >= 0.0_rp ) then
                 kfl_fixno_lev(1,ipoin) = -8              ! Outflow
              else
                 kfl_fixno_lev(1,ipoin) =  8              ! Inflow
              end if
           end if
        end if
     end do
     if ( kfl_conbc_lev == 0 ) then
        print*,'----------updbcs-----------------'
        relax = 0.5_rp
           do ipoin = 1,npoin
              if( kfl_funno_lev(ipoin) < 0 ) then 
                 ifunc = -kfl_funno_lev(ipoin)            
                 call ker_space_time_function(&
                      ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,lefun)
                 print*,'ipoin,fleve',ipoin,fleve(ipoin,1),kfl_fixno_lev(1,ipoin)
                 !if(kfl_fixno_lev(1,ipoin)==11)fleve(ipoin,1) = lefun
                 if(kfl_fixno_lev(1,ipoin)==11)fleve(ipoin,1) = fleve(ipoin,1) + lefun
                 if(kfl_fixno_lev(1,ipoin)==12)fleve(ipoin,1) = fleve(ipoin,1) + lefun * (1.0_rp - relax)
                 print*,'ipoin,fleve',ipoin,fleve(ipoin,1),lefun
                 !if( kfl_timei_lev /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
                 !   do idime = 1,ndime
                 !      veold = veloc(idime,ipoin,ncomp_nsi)
                 !      bvess_nsi(idime,ipoin,1) = 0.50_rp*(vefun(idime)+veold)
                 !   end do
                 !else
                 !   do idime = 1,ndime
                 !      bvess_nsi(idime,ipoin,1) = vefun(idime)
                 !   end do
                 !end if
              end if
           end do
        end if
     end if

end subroutine lev_updbcs

