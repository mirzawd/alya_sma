!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_iniscv.f90
!> @author  Jazmin Aguado-Sierra
!> @brief   Initial condition setup for Paci 2013 Ventricular-like model
!> @details provides unique initial conditions to the FEM
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_iniscv(ipoin)
!subroutine exm_iniscv(kmodel_ipoin,ipoin)

  use      def_master
  use      def_domain
  use      def_elmtyp
  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: ipoin
  !integer(ip)   :: mat
  !ituss = 0_ip  !!%endo = 1, epi = 0, M = 2  

  !if(INOTMASTER) then
  ! do imate= 1,nmate_exm
  !  clmat=kfl_cellmod(imate)
  !  if(clmat == 7) then
  !     n=imate
  !  end if        
  ! end do
      
  !   mat=kmodel_ipoin-3    
  
  elmag(ipoin,1:3) = -70.0_rp     ! resting voltage
  vicel_exm(1:15,ipoin) = 0.0_rp  !all currents
  vauxi_exm(1,ipoin,1:3) = 0.75_rp  !h0    
  vauxi_exm(2,ipoin,1:3) =  0.75_rp  !j0
  vauxi_exm(3,ipoin,1:3) = 0.0_rp   !m0
  vauxi_exm(4,ipoin,1:3) = 0.0_rp   !d0    
  vauxi_exm(5,ipoin,1:3) = 1.0_rp   !f_Ca0
  vauxi_exm(6,ipoin,1:3) = 1.0_rp   !f10
  vauxi_exm(7,ipoin,1:3) = 1.0_rp   !f20  
  vauxi_exm(8,ipoin,1:3) = 0.0_rp  !r0
  vauxi_exm(9,ipoin,1:3) = 1.0_rp   !q0
  vauxi_exm(10,ipoin,1:3) = 0.0_rp    !xr10  
  vauxi_exm(11,ipoin,1:3) = 1.0_rp  !xr20
  vauxi_exm(12,ipoin,1:3) = 0.0_rp    !Xs0
  vauxi_exm(13,ipoin,1:3) = 0.1_rp    !Xf0
  vauxi_exm(14,ipoin,1:3) = 1.0_rp    !g0
  vconc(1,ipoin,1:3) = 0.0002_rp   !Ca_i0
  vconc(2,ipoin,1:3) = 0.3_rp    !Ca_SR0
  vconc(3,ipoin,1:3) = 10.0_rp    !Na_i0
end subroutine exm_iniscv
