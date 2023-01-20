!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup SolidzMatrixAssembly
!> @ingroup    Solidz
!> @{
!> @file    sld_minmax.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute some min/max/sum values
!> @details Compute some valeus ot output in the convergence file
!>
!> @} 
!------------------------------------------------------------------------

subroutine sld_minmax(cpu_ass_ave,cpu_ass_max,cpu_sol_ave,cpu_sol_max,displ_max)

  use def_kintyp,         only : ip,rp
  use def_master,         only : unkno
  use def_domain,         only : npoin,ndime
  use def_solidz,         only : cpu_ass_sol_sld
  use def_solidz,         only : ndofn_sld, volum_sld
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE
  use mod_communications, only : PAR_SUM
  implicit none

  real(rp),    intent(out) :: cpu_ass_ave
  real(rp),    intent(out) :: cpu_ass_max
  real(rp),    intent(out) :: cpu_sol_ave
  real(rp),    intent(out) :: cpu_sol_max
  real(rp),    intent(out) :: displ_max
  integer(ip)              :: ipoin,idofn
  real(rp)                 :: rmaxi(5),raver(5)
  !
  ! DISPL: Max displacement
  !
  displ_max = -huge(1.0_rp)
  do ipoin = 1,npoin
     idofn = ( ipoin-1 ) * ndofn_sld
     displ_max = max(displ_max,dot_product(unkno(idofn+1:idofn+ndime),unkno(idofn+1:idofn+ndime)))
  end do
  displ_max = sqrt(max(0.0_rp,displ_max))
  !
  ! Call to parall service
  !
  rmaxi(1)    = displ_max
  rmaxi(2)    = cpu_ass_sol_sld(1)
  rmaxi(3)    = cpu_ass_sol_sld(2)
  call PAR_MAX(3_ip,rmaxi)
  displ_max   = rmaxi(1)
  cpu_ass_max = rmaxi(2)
  cpu_sol_max = rmaxi(3)
  
  raver(1)    = cpu_ass_sol_sld(1)
  raver(2)    = cpu_ass_sol_sld(2)
  call PAR_AVERAGE(2_ip,raver)
  cpu_ass_ave = raver(1)
  cpu_sol_ave = raver(2)

  call PAR_SUM(volum_sld(1),'IN MY CODE')
  call PAR_SUM(volum_sld(2),'IN MY CODE')

end subroutine sld_minmax

