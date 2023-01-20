!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine proaxy(nbvar,a,xx,yy,sumxx)
  !------------------------------------------------------------------------
  !****f* solite/proaxy
  ! NAME 
  !    proaxy
  ! DESCRIPTION
  !    This routine computes the sclara product of two vectors:
  !    SUMXX=XX^t.YY
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr
  use def_domain, only     :  npoin
  implicit none
  integer(ip), intent(in)  :: nbvar
  real(rp),    intent(in)  :: a,xx(*),yy(*)
  real(rp),    intent(out) :: sumxx
  real(rp),    target      :: dummr_par(2) 
  integer(ip)              :: ii
  real(rp)                 :: numer,denom,dummr

  numer = 0.0_rp
  denom = 0.0_rp

  if(kfl_paral==-1) then
     !
     ! Sequential
     !
     do ii=1,nbvar*npoin
        dummr = xx(ii) - a * yy(ii)
        numer = numer + dummr  * dummr
        denom = denom + xx(ii) * xx(ii)
     end do

  else if(kfl_paral>=1) then
     !
     ! Parallel: Slaves
     !
     if(nbvar==1) then
        do ii=1,npoi1
           dummr = xx(ii) - a * yy(ii)
           numer = numer + dummr  * dummr
           denom = denom + xx(ii) * xx(ii)
        end do
        do ii=npoi2,npoi3
           dummr = xx(ii) - a * yy(ii)
           numer = numer + dummr  * dummr
           denom = denom + xx(ii) * xx(ii)
        end do
     else
        do ii=1,nbvar*npoi1
           dummr = xx(ii) - a * yy(ii)
           numer = numer + dummr  * dummr
           denom = denom + xx(ii) * xx(ii)
        end do
        do ii=(npoi2-1)*nbvar+1,npoi3*nbvar
           dummr = xx(ii) - a * yy(ii)
           numer = numer + dummr  * dummr
           denom = denom + xx(ii) * xx(ii)
        end do
     end if
  end if

  if(kfl_paral>=0) then
     !
     ! Parallel: reduce sum
     !
     nparr        =  2
     dummr_par(1) =  numer
     dummr_par(2) =  denom
     parre        => dummr_par
     call par_operat(3_ip)           
     numer        =  dummr_par(1)
     denom        =  dummr_par(2)
  end if

  sumxx = sqrt(numer/denom)

end subroutine proaxy

