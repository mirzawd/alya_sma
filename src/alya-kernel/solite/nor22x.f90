!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nor22x(nbvar,nbnodes,xx,sumxx)
  !------------------------------------------------------------------------
  !****f* solite/nor22x
  ! NAME 
  !    prodxy
  ! DESCRIPTION
  !    This routine computes the square of the 2-norm of a vector XX:
  !    SUMXX = ||XX||^2_2 = ( sum_i |XX_i|^2 )
  ! INPUT
  !    NBNODES .................. Number of nodes
  !    NBVAR .................... Number of variables per node
  !    XX(NBNODES*NBVAR) ........ Vector
  ! OUTPUT
  !    SUMXX .................... Eucledian norm
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr,npari,nparc
  implicit none
  integer(ip), intent(in)  :: nbvar,nbnodes
  real(rp),    intent(in)  :: xx(*)
  real(rp),    intent(out) :: sumxx
  real(rp),    target      :: dummr_par(1) 
  integer(ip)              :: kk,ii

  sumxx = 0.0_rp

  if(kfl_paral==-1) then
     !
     ! Sequential
     !
     do kk=1,nbvar*nbnodes
        sumxx = sumxx + xx(kk) * xx(kk)
     end do

  else if(kfl_paral>=1) then
     !
     ! Parallel: slaves
     !
     if(nbvar==1) then
        do ii=1,npoi1
           sumxx  = sumxx  + xx(ii) * xx(ii)
        end do
        do ii=npoi2,npoi3
           sumxx  = sumxx  + xx(ii) * xx(ii)
        end do

     else
        do ii=1,npoi1*nbvar
           sumxx  = sumxx  + xx(ii) * xx(ii)
        end do
        do ii=(npoi2-1)*nbvar+1,npoi3*nbvar
           sumxx  = sumxx  + xx(ii) * xx(ii)
        end do
     end if
  end if
  !
  ! Parallel: reduce sum
  !
  npari=0
  nparc=0
  if(kfl_paral>=0) then
     nparr        =  1
     dummr_par(1) =  sumxx
     parre        => dummr_par
     call par_operat(3_ip)
     sumxx        =  dummr_par(1)
  end if

end subroutine nor22x
