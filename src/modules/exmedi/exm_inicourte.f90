!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exmedi_inicourte.f90
!> @author  Eva Casoni 
!> @date    10/12/2020
!> @brief   Initial condition setup for Courtemanche atrial  model
!> @details 
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_inicourte(ipoin,mat)

  use def_master
  use def_domain
  use def_elmtyp
  use def_exmedi
  use mod_exm_cellmodel


  implicit none
  integer(ip), intent(in) :: ipoin, mat
  integer(ip)             :: ituss_exm
  !integer(ip) ::  mat

  !
  ! Load initial conditions for the potentials
  !
  !mat=nodemat(ituss_exm)
  ituss_exm = int(celty_exm(ipoin), kind=ip)

  elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)              !Membrance (potential), V

  vauxi_exm(1:15,ipoin,1) = vauxi_exm_initial(1:15,ituss_exm,mat)          !Gate variables
  vauxi_exm(1:15,ipoin,2) = vauxi_exm_initial(1:15,ituss_exm,mat)          
  vauxi_exm(1:15,ipoin,3) = vauxi_exm_initial(1:15,ituss_exm,mat)          
  
  vconc(1,ipoin,1) = vconc_initial(1,ituss_exm,mat)              !Variable Cai
  vconc(2,ipoin,1) = vconc_initial(2,ituss_exm,mat)              !Variable Carel
  vconc(3,ipoin,1) = vconc_initial(3,ituss_exm,mat)              !Variable Caup
  vconc(4,ipoin,1) = vconc_initial(4,ituss_exm,mat)              !Variable Ki
  vconc(5,ipoin,1) = vconc_initial(5,ituss_exm,mat)              !Variable Nai
  vconc(1,ipoin,2) = vconc_initial(1,ituss_exm,mat)              !Variable Cai
  vconc(2,ipoin,2) = vconc_initial(2,ituss_exm,mat)              !Variable Carel
  vconc(3,ipoin,2) = vconc_initial(3,ituss_exm,mat)              !Variable Caup
  vconc(4,ipoin,2) = vconc_initial(4,ituss_exm,mat)              !Variable Ki
  vconc(5,ipoin,2) = vconc_initial(5,ituss_exm,mat)              !Variable Nai
  vconc(1,ipoin,3) = vconc_initial(1,ituss_exm,mat)              !Variable Cai
  vconc(2,ipoin,3) = vconc_initial(2,ituss_exm,mat)              !Variable Carel
  vconc(3,ipoin,3) = vconc_initial(3,ituss_exm,mat)              !Variable Caup
  vconc(4,ipoin,3) = vconc_initial(4,ituss_exm,mat)              !Variable Kio
  vconc(5,ipoin,3) = vconc_initial(5,ituss_exm,mat)              !Variable Nai

  vicel_exm(1:18,ipoin) = 0.0_rp



end subroutine exm_inicourte

