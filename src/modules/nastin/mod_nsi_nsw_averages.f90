!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_nsw_averages.f90
!> @author  Herbert Owen
!> @brief   Obtains some values and averages needed to apply wall law increasing viscosity in the first element and using no slip 
!> @details - obtains: 
!>          - 
!> @} 
!-----------------------------------------------------------------------
module mod_nsi_nsw_averages
  use def_kintyp,                   only : ip,rp
  use def_nastin,                   only : massb_nsi,notra_nsi,vafor_nsi,avvaf_nsi,avntr_nsi,avgtr_nsi
  use def_parame,                   only : zero
  use def_master,                   only : kfl_paral

  implicit none
  private
  public :: nsi_nsw_averages

contains
  subroutine nsi_nsw_averages(avwei)
    use def_master,                   only : INOTMASTER
    use def_domain,                   only : ndime,npoin,lpoty
    use def_master,                   only : gevec

    implicit none
    real(rp),intent(in)               :: avwei

    integer(ip)   :: ipoin,idime,ibopo
    real(rp)      :: auvar,augra

    if( INOTMASTER ) then     ! pensar que pasa con el master
       !
       ! time average tangential force calculated with teh velocidty gradientes
       ! it is very similar to AVTAN but it is a running averge instead of a real average
       !
       call memgen(zero,ndime,npoin)
       call nsi_memphy(4_ip)
       call nsi_denvis()
       call nsi_outtan(0_ip)    ! here I want gradient based traction - It comes in gevec  
       call nsi_memphy(-4_ip)
       !
       ! running Time averages 
       !
       do ipoin = 1,npoin
          ibopo=lpoty(ipoin)
          if(ibopo > 0_ip) then 
             auvar = 0.0_rp
             augra = 0.0_rp
             do idime=1,ndime
                if (abs(massb_nsi(ipoin)) > 1.0e-30_rp) notra_nsi(idime,ipoin) =  &
                     - vafor_nsi(idime,ipoin) / massb_nsi(ipoin) ! added negative sign so that it has the same dir as tange
                avvaf_nsi(idime,ipoin) = (1.0_rp-avwei) * avvaf_nsi(idime,ipoin)  &
                     + avwei  * vafor_nsi(idime,ipoin)
                avntr_nsi(idime,ipoin) = (1.0_rp-avwei) * avntr_nsi(idime,ipoin) &
                     + avwei  * notra_nsi(idime,ipoin)  ! variational based - reaction al nodes
                avgtr_nsi(idime,ipoin) = (1.0_rp-avwei) * avgtr_nsi(idime,ipoin) &
                     + avwei  * gevec(idime,ipoin)      ! calculated from velocity gradients
                auvar = auvar + avntr_nsi(idime,ipoin) * avntr_nsi(idime,ipoin) 
                augra = augra + avgtr_nsi(idime,ipoin) * avgtr_nsi(idime,ipoin)
             end do
             auvar = sqrt(auvar)
             augra = sqrt(augra)
          else
             avvaf_nsi(:,ipoin) = 0.0_rp
             avntr_nsi(:,ipoin) = 0.0_rp
             avgtr_nsi(:,ipoin) = 0.0_rp
          end if
       end do
       call memgen(2_ip,ndime,npoin)

    end if

  end subroutine nsi_nsw_averages

end module mod_nsi_nsw_averages
