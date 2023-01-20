!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_exm_ecg
   use def_kintyp_basic, only : i1p, ip, rp, lg 
   implicit none

   private

   public :: exm_ecg_exchange                  ! Initial exchange of related variables, to put into sendat. Call once!
   public :: exm_ecg_allocate                  ! Allocate memory for stimulus table
contains

subroutine exm_ecg_exchange()
   use mod_exchange, only : exchange_add, exchange_end, exchange_init
   use def_exmedi, only : pseudecg_exm, ecg_points, nrootecg_exm
   implicit none

   integer(ip) :: i


   if (nrootecg_exm>0_ip) then
      call exchange_init()
      call exchange_add( pseudecg_exm )
      
      do i=1,size(ecg_points, kind=ip)
         call exchange_add( ecg_points(i) % coords )
         call exchange_add( ecg_points(i) % label  )
      end do

      call exchange_end()
   end if
end subroutine exm_ecg_exchange

subroutine exm_ecg_allocate()
   use def_exmedi,     only : pseudecg_exm, ecg_points, nrootecg_exm
   use mod_exm_memory, only : memory_alloca
   use def_master,     only : mem_modul, modul
   implicit none

   if (nrootecg_exm>0_ip) then
      call memory_alloca(mem_modul(1:2,modul),'pseudecg_exm', 'mod_exm_ecg', pseudecg_exm, nrootecg_exm)
      pseudecg_exm = 0.0_rp

      call memory_alloca(mem_modul(1:2,modul),'ecg_points', 'mod_exm_ecg', ecg_points, nrootecg_exm)
   end if
end subroutine exm_ecg_allocate



end module mod_exm_ecg
