!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine prodxy(nbvar,nbnodes,xx,yy,sumxy)
  !------------------------------------------------------------------------
  !****f* solite/prodxy
  ! NAME 
  !    prodxy
  ! DESCRIPTION
  !    This routine computes the sclara product of two vectors:
  !    SUMXY=XX^t.YY
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only : ip,rp
  use def_master,         only : kfl_paral,npoi3
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip), intent(in)  :: nbvar,nbnodes
  real(rp),    intent(in)  :: xx(*),yy(*)
  real(rp),    intent(out) :: sumxy
#ifdef BLAS
  real(rp) :: DDOT
  external DDOT
#endif

  if( kfl_paral == -1 ) then
#ifdef BLAS
     sumxy = DDOT(nbvar*nbnodes,xx,1_ip,yy,1_ip)
#else
     sumxy = dot_product(xx(1:nbvar*nbnodes),yy(1:nbvar*nbnodes))
#endif
     
  else if( kfl_paral >= 1 ) then
#ifdef BLAS
     sumxy = DDOT(nbvar*npoi3,xx,1_ip,yy,1_ip)
#else  
     sumxy = dot_product(xx(1:nbvar*npoi3),yy(1:nbvar*npoi3))
#endif
  end if
  
  if( kfl_paral >= 0 ) call PAR_SUM(sumxy)

end subroutine prodxy
