!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine norm1x(nbvar,xx,sumxx)
  !------------------------------------------------------------------------
  !****f* solite/norm2x
  ! NAME 
  !    prodxy
  ! DESCRIPTION
  !    This routine computes the 1-norm of a vector XX:
  !    SUMXX = ||XX||_1 = ( sum_i |XX_i| )
  ! INPUT
  !    NBVAR .................... Number of variables per node
  !    XX(NBVAR,NPOIN) .......... Vector
  ! OUTPUT
  !    SUMXX .................... 1-norm
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  npoin
  use def_master, only     :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr
  implicit none
  integer(ip), intent(in)  :: nbvar
  real(rp),    intent(in)  :: xx(*)
  real(rp),    intent(out) :: sumxx
  real(rp),    target      :: dummr_par(1) 
  integer(ip)              :: kk,ii

  sumxx = 0.0_rp

  if(kfl_paral==-1) then
     !
     ! Sequential
     !
     do kk=1,nbvar*npoin
        sumxx = sumxx + abs(xx(kk))
     end do

  else if(kfl_paral>=1) then
     !
     ! Parallel: Slaves
     !
     if(nbvar==1) then
        do ii=1,npoi1
           sumxx  = sumxx  + abs(xx(ii))
        end do
        do ii=npoi2,npoi3
           sumxx  = sumxx  + abs(xx(ii))
        end do
     else 
        do ii=1,npoi1*nbvar
           sumxx  = sumxx  + abs(xx(ii))
        end do
        do ii=(npoi2-1)*nbvar+1,npoi3*nbvar
           sumxx  = sumxx  + abs(xx(ii))
        end do
     end if
  end if

  if(kfl_paral>=0) then
     !
     ! Master and Slaves
     !
     nparr        =  1
     dummr_par(1) =  sumxx
     parre        => dummr_par
     call par_operat(3_ip)           
     sumxx        =  dummr_par(1)
  end if

end subroutine norm1x
