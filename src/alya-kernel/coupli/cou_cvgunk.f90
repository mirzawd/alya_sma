!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling arrays
!> @file    def_coupli.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Output coupling convergence
!> @details Output coupling convergence
!> @}
!------------------------------------------------------------------------

subroutine cou_cvgunk()
  use def_parame
  use def_master
  use def_domain
  use def_coupli
  use mod_couplings, only : THERE_EXISTS_A_ZONE_COUPLING
  use mod_iofile,    only : iofile_flush_unit
  implicit none
  integer(ip), save :: ipass = 0
  integer(ip)       :: icoup
  integer(4)        :: lun_coupl_cvg4
  real(rp)          :: retot
  !
  ! Write convergence
  !
  if( INOTSLAVE .and. mcoup > 0 .and. THERE_EXISTS_A_ZONE_COUPLING() ) then 
     lun_coupl_cvg4 = int(lun_coupl_cvg,4)
     retot = 0.0_rp 
     do icoup = 1,mcoup
        retot = retot + coupling_type(icoup) % resid(1)
     end do
     if( ipass == 0 ) write(lun_coupl_cvg4,100)
     retot = retot / real(mcoup,rp)
     write(lun_coupl_cvg4,101) &
          ittim,iblok,cutim,retot,(coupling_type(icoup) % resid(1),icoup=1,mcoup)
     call iofile_flush_unit(lun_coupl_cvg4)
  end if

  ipass = 1     

  !----------------------------------------------------------------------
  !
  ! Formats
  !
  !----------------------------------------------------------------------

100 format('# --| ALYA convergence  ' ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --|  1. Time Step             2. Block number           3. Current time      ',/,&
       & '# --|  4. Global residual       5. Resid. coupling 1      6. Resid coupling 2  ',/,&
       & '# --|  7. Resid coupling 3      8. => All other residuals                      ',/,&
       & '# ','             1','             2','             3','             4'&
       &     ,'             5','             6','             7','             8') 
101 format(4x,i12,2x,i12,42(2x,e16.8e3))

end subroutine cou_cvgunk
