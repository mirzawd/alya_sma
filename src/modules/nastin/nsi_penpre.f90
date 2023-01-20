!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_penpre.f90
!> @author  Guillaume Houzeaux
!> @brief   Pressure penalization
!> @details Penalize the pressure equation as:
!>          \verbatim
!>          App.p^i+1 + eps*diag(App).p^i + Apu.u = bp + eps*diag(App) p^{i-1}
!>          Q <= eps*diag(App)
!>          \endverbatim
!>          where eps is a user-defined constant
!> @} 
!------------------------------------------------------------------------

subroutine nsi_penpre(App,Q,bp)
  use def_kintyp
  use def_domain
  use def_master
  use def_nastin

  implicit none
  real(rp),    intent(out)   :: App(nzdom)
  real(rp),    intent(out)   :: Q(nzdom)
  real(rp),    intent(out)   :: bp(npoin)
  integer(ip)                :: ipoin,jpoin,izdom
  real(rp)                   :: xdiag
  
  do ipoin = 1,npoin

     izdom = r_dom(ipoin) - 1
     jpoin = 0
     do while( jpoin /= ipoin )
        izdom = izdom + 1
        jpoin = c_dom(izdom)
     end do 

     xdiag      = App(izdom) * prepe_nsi
     App(izdom) = App(izdom) + xdiag
     Q(izdom)   = Q(izdom)   + xdiag
     bp(ipoin)  = bp(ipoin)  + xdiag * press(ipoin,1)

  end do

end subroutine nsi_penpre

