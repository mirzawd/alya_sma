!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine norm2x(nbvar,xx,sumxx)
  !------------------------------------------------------------------------
  !****f* solite/norm2x
  ! NAME 
  !    prodxy
  ! DESCRIPTION
  !    This routine computes the 2-norm (Eucledian norm) of a vector XX:
  !    SUMXX = ||XX||_2 = ( sum_i |XX_i|^2 )^{1/2}
  ! INPUT
  !    NBVAR .................... Number of variables per node
  !    XX(NBVAR,NPOIN) .......... Vector
  ! OUTPUT
  !    SUMXX .................... Eucledian norm
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only : ip,rp
  use def_domain,         only : npoin
  use def_master,         only : kfl_paral,npoi3
!  use def_master,         only : npoi1,npoi2
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip), intent(in)  :: nbvar
  real(rp),    intent(in)  :: xx(*)
  real(rp),    intent(out) :: sumxx
!  integer(ip)              :: ii,jj
#ifdef BLAS2
  real(rp) :: DDOT
  external DDOT
#endif

  if( kfl_paral == -1 ) then
     !
     ! Sequential
     !
#ifdef BLAS2
     sumxx = DDOT(nbvar*npoin,xx,1_ip,xx,1_ip)
#else
     sumxx = dot_product(xx(1:nbvar*npoin),xx(1:nbvar*npoin))
#endif

  else if( kfl_paral >= 1 ) then
     !
     ! Parallel: slaves
     !
#ifdef BLAS2
     sumxx = DDOT(nbvar*npoi3,xx,1_ip,xx,1_ip)
#else
     sumxx = dot_product(xx(1:nbvar*npoi3),xx(1:nbvar*npoi3))
     !sumxx = dot_product(xx(1:nbvar*npoi1),xx(1:nbvar*npoi1))
#endif
     !ii    = (npoi2-1)*nbvar+1
     !jj    = npoi3*nbvar
     !sumxx = sumxx + dot_product(xx(ii:jj),xx(ii:jj))
  end if

  if( kfl_paral >= 0 ) call PAR_SUM(sumxx)
  sumxx = sqrt(sumxx)

end subroutine norm2x
