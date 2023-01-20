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
!> @brief   Open coupling files
!> @details Open coupling files
!> @}
!------------------------------------------------------------------------

subroutine cou_openfi()
  use def_parame
  use def_master
  use def_domain
  use def_coupli
  use mod_iofile
  implicit none
  character(150) :: fil_coupl_res,fil_coupl_dat,fil_coupl_cvg

  if( INOTSLAVE ) then

     fil_coupl_dat = adjustl(trim(namda)) // '.cou.dat'
     fil_coupl_res = adjustl(trim(namda)) // '.cou.log'
     fil_coupl_cvg = adjustl(trim(namda)) // '.cou.cvg'

     call iofile(7_ip,lun_coupl_dat,trim(fil_coupl_dat),'COUPLING DATA','old')
     if( .not. file_opened ) lun_coupl_dat = 0

     if( lun_coupl_dat /= 0 ) then
        if( kfl_rstar == 2 ) then
           call iofile(0_ip,lun_coupl_res,fil_coupl_res,'COUPLING RESULTS'    ,'old','formatted','append')
           call iofile(0_ip,lun_coupl_cvg,fil_coupl_cvg,'COUPLING CONVERGENCE','old','formatted','append')
        else
           call iofile(0_ip,lun_coupl_res,fil_coupl_res,'COUPLING RESULTS')
           call iofile(0_ip,lun_coupl_cvg,fil_coupl_cvg,'COUPLING CONVERGENCE')
        end if
     end if

  end if

end subroutine cou_openfi
