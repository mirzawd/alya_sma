!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outres()

!------------------------------------------------------------------------
!****f* Nastin/nsi_outres
! NAME 
!    nsi_outres
! DESCRIPTION
!    This routine output the velocity residual
! USES
! USED BY
!    nsi_output
!------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin
  use mod_communications, only : PAR_MAX
  implicit none
  integer(ip) :: ipoin,idime
  real(rp)    :: rnorm,resid(2),rmax
  !
  ! Compute residual
  !
  rnorm = 0.0_rp
  
  do ipoin = 1,npoin
     resid = 0.0_rp
     do idime = 1,ndime
        resid(1) = resid(1) + (veloc(idime,ipoin,1)-veold_nsi(idime,ipoin))**2
        resid(2) = resid(2) + veloc(idime,ipoin,1)**2
     end do
     resid        = sqrt(resid)
     rnorm        = max(rnorm,resid(2))
     gesca(ipoin) = resid(1)
  end do

  call PAR_MAX(rnorm)

  rmax = 0.0_rp
  do ipoin = 1,npoin
     gesca(ipoin) = gesca(ipoin)/(rnorm+zeror)
     rmax = max(rmax,gesca(ipoin))
  end do

end subroutine nsi_outres
