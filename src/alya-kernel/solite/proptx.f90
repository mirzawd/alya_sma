!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!> proptx.f90
!> @file proptx.f90 
!> @fn proptx 
!> Text automatically added for recognition of doxygen parser
!>
subroutine proptx(nbnodes,nbvar,xx,yy,sumxx)

  !------------------------------------------------------------------------
  ! Sources/kernel/solite/proptx.f90
  ! NAME 
  !    proptx
  ! DESCRIPTION
  !    This routine computes the sclara product of two complex vectors:
  !                      SUMXX = <XX,YY>
  ! INPUT ARGUMENTS
  !    NBNODES .................. Number of nodes
  !    NBVAR .................... Number of variables per node
  !    XX(NBVAR*NBNODES) ........ Complex vector #1
  !    YY(NBVAR*NBNODES) ........ Complex vector #2
  ! OUTPUT ARGUMENTS
  !    SUMXX .................... Scalar (dot) product 
  ! USES
  ! USED BY 
  !------------------------------------------------------------------------


  use def_kintyp, only : ip,rp
  use def_master, only : kfl_paral,npoi1,npoi2,npoi3
  use mod_communications, only : PAR_SUM
! Declaration statements
  implicit none

! Dummy arguments
  integer(ip), intent(in)  :: nbnodes,nbvar
  complex(rp), intent(in)  :: xx(*),yy(*)
  complex(rp), intent(out) :: sumxx

! Local variables
  integer(ip)              :: kk,temp

!--------------------------------------------------------------------------------------------------------
! Product: sumxx = x^T * y 
!--------------------------------------------------------------------------------------------------------

  sumxx = (0.0_rp,0.0_rp)
 
  if ( kfl_paral == -1_ip ) then
  !Sequential
	do kk = 1,nbvar*nbnodes

		sumxx = sumxx + yy(kk) * xx(kk) 
		
	end do

  else if ( kfl_paral >= 1_ip ) then
  !Parallel: slaves
  	do kk = 1,nbvar*npoi1

        	sumxx = sumxx + yy(kk) * xx(kk) 
		
  	end do     
        temp = nbvar*(npoi2 - 1_ip) + 1_ip
  	do kk = temp,nbvar*npoi3

        	sumxx = sumxx + yy(kk) * xx(kk) 
		
	end do

  end if

  call PAR_SUM(sumxx)
 
end subroutine proptx


