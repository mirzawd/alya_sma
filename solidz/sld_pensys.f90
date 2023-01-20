!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_penpre.f90
!> @author  Mariano Vazquez
!> @brief   Correct the system matrix to include penalization by the diagonal 
!> @details Penalize the pressure equation as:
!>          \verbatim
!>          [A + eps*diag(A)] DU = b
!>          \endverbatim
!>          where eps is a user-defined constant
!> @} 
!------------------------------------------------------------------------
subroutine sld_pensys(amatrix)
  use def_kintyp
  use def_domain
  use def_master
  use def_solidz

  implicit none
  real(rp)                   :: amatrix(nzdom)
  integer(ip)                :: ipoin,jpoin,izdom
  real(rp)                   :: xdiag
  
  do ipoin = 1,npoin

     izdom = r_dom(ipoin) - 1
     jpoin = 0
     do while( jpoin /= ipoin )
        izdom = izdom + 1
        jpoin = c_dom(izdom)
     end do 

     xdiag      = amatrix(izdom) * factor_penal_sld
     amatrix(izdom) = amatrix(izdom) + xdiag

  end do

end subroutine sld_pensys
