!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine prxyxz(nbvar,xx,yy,zz,sumxy,sumxz)
  !------------------------------------------------------------------------
  !****f* solite/prxyxz
  ! NAME 
  !    prxyxz
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
  real(rp),    intent(in)  :: xx(*),yy(*),zz(*)
  real(rp),    intent(out) :: sumxy,sumxz
  real(rp),    target      :: dummr_par(2) 
  integer(ip)              :: kk,ii

  sumxy = 0.0_rp
  sumxz = 0.0_rp

  if(kfl_paral==-1) then
     !
     ! Sequential
     !
     do kk=1,nbvar*npoin
        sumxy = sumxy + xx(kk) * yy(kk)
        sumxz = sumxz + xx(kk) * zz(kk)
     end do

  else if(kfl_paral>=1) then
     !
     ! Parallel: Slaves
     !
     if(nbvar==1) then
        do ii=1,npoi1
           sumxy  = sumxy  + xx(ii) * yy(ii)
           sumxz  = sumxz  + xx(ii) * zz(ii)
        end do
        do ii=npoi2,npoi3
           sumxy  = sumxy  + xx(ii) * yy(ii)
           sumxz  = sumxz  + xx(ii) * zz(ii)
        end do
     else
        do ii=1,nbvar*npoi1
           sumxy = sumxy + xx(ii) * yy(ii)
           sumxz = sumxz + xx(ii) * zz(ii)
        end do
        do ii=(npoi2-1)*nbvar+1,npoi3*nbvar
           sumxy  = sumxy  + xx(ii) * yy(ii)
           sumxz  = sumxz  + xx(ii) * zz(ii)
        end do
     end if
  end if

  if(kfl_paral>=0) then
     !
     ! Parallel: reduce sum
     !
     nparr        =  2
     dummr_par(1) =  sumxy
     dummr_par(2) =  sumxz
     parre        => dummr_par
     call par_operat(3_ip)           
     sumxy        =  dummr_par(1)
     sumxz        =  dummr_par(2)
  end if

end subroutine prxyxz

 
