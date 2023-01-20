!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!>
!> @addtogroup Nastin
!> @{
!> @file    nsi_minmax.f90
!> @author  houzeaux
!> @date    2018-08-27
!> @brief   Statistics
!> @details Compute some statistic on timings, min and max values, etc.
!> @}
!> 
!-----------------------------------------------------------------------

subroutine nsi_minmax(&
     cpu_ass_ave,cpu_ass_max,cpu_bou_ave,cpu_bou_max,&
     cpu_sol_ave,cpu_sol_max,cpu_sgs_ave,cpu_sgs_max)
  
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE
  implicit none
  
  real(rp),    intent(out) :: cpu_ass_ave
  real(rp),    intent(out) :: cpu_ass_max
  real(rp),    intent(out) :: cpu_bou_ave
  real(rp),    intent(out) :: cpu_bou_max
  real(rp),    intent(out) :: cpu_sol_ave
  real(rp),    intent(out) :: cpu_sol_max
  real(rp),    intent(out) :: cpu_sgs_ave
  real(rp),    intent(out) :: cpu_sgs_max
  integer(ip), parameter   :: nmima=3
  integer(ip)              :: ipoin,idime,ndofn,idofn
  integer(ip)              :: ifac1,ifac2
  real(rp)                 :: rmini(5),rmaxi(7),raver(4)
  real(rp)                 :: vnorm

  vemin_nsi =  huge(1.0_rp)
  vemax_nsi = -huge(1.0_rp)
  prmin_nsi =  huge(1.0_rp)
  prmax_nsi = -huge(1.0_rp)
  if( NSI_MONOLITHIC ) then
     ndofn = ndime+1
     ifac1 = ndofn
     ifac2 = 0
  else
     ndofn = ndime
     ifac1 = 1
     ifac2 = npoin * ndofn
  end if
  !
  ! VEMIN and VEMAX: Min and max pressure
  ! PRMIN and PRMAX: Min and max pressure
  !
  do ipoin = 1,npoin
     idofn = ( ipoin-1 ) * ndofn
     vnorm = 0.0_rp
     do idime = 1,ndime
        idofn = idofn + 1
        vnorm = vnorm + unkno(idofn) * unkno(idofn)
     end do
     idofn = ipoin * ifac1 + ifac2  
     vemax_nsi = max(vemax_nsi,vnorm)
     vemin_nsi = min(vemin_nsi,vnorm)
     prmax_nsi = max(prmax_nsi,unkno(idofn))
     prmin_nsi = min(prmin_nsi,unkno(idofn))
  end do
  vemin_nsi = sqrt(max(0.0_rp,vemin_nsi))
  vemax_nsi = sqrt(max(0.0_rp,vemax_nsi))
  !
  ! Call to parall service
  !
  rmini(1)  =  tamin_nsi
  rmini(2)  =  vemin_nsi
  rmini(3)  =  prmin_nsi
  call PAR_MIN(3_ip,rmini,'IN MY CODE')
  tamin_nsi =  rmini(1)
  vemin_nsi =  rmini(2)
  prmin_nsi =  rmini(3)

  rmaxi(1)    = tamax_nsi
  rmaxi(2)    = vemax_nsi
  rmaxi(3)    = prmax_nsi
  rmaxi(4)    = cpu_ass_sol_nsi(1)
  rmaxi(5)    = cpu_ass_sol_nsi(2)
  rmaxi(6)    = cpu_ass_sol_nsi(3)
  rmaxi(7)    = cpu_ass_sol_nsi(4)
  call PAR_MAX(7_ip,rmaxi,'IN MY CODE')
  tamax_nsi   = rmaxi(1)
  vemax_nsi   = rmaxi(2)
  prmax_nsi   = rmaxi(3)
  cpu_ass_max = rmaxi(4)
  cpu_bou_max = rmaxi(7)
  cpu_sol_max = rmaxi(5)
  cpu_sgs_max = rmaxi(6)

  raver(1)    = cpu_ass_sol_nsi(1)
  raver(2)    = cpu_ass_sol_nsi(2)
  raver(3)    = cpu_ass_sol_nsi(3)
  raver(4)    = cpu_ass_sol_nsi(4)
  call PAR_AVERAGE(4_ip,raver,'IN MY CODE')
  cpu_ass_ave = raver(1)
  cpu_bou_ave = raver(4)
  cpu_sol_ave = raver(2)
  cpu_sgs_ave = raver(3)

end subroutine nsi_minmax

