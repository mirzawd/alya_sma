!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine proxyz(nbvar,nbnodes,xx,yy,zz,sumxx)
  !------------------------------------------------------------------------
  !****f* solite/proxyz
  ! NAME 
  !    proxyz
  ! DESCRIPTION
  !    This routine computes the sclara product of two vectors:
  !    SUMXX=XX^t.YY
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr
  implicit none
  integer(ip), intent(in)  :: nbvar,nbnodes
  real(rp),    intent(in)  :: xx(*),yy(*),zz(*)
  real(rp),    intent(out) :: sumxx
  real(rp),    target      :: dummr_par(2) 
  integer(ip)              :: kk,ii

  sumxx = 0.0_rp

  if( kfl_paral == -1 ) then
     !
     ! Sequential
     !
     do kk=1,nbvar*nbnodes
        sumxx = sumxx + xx(kk) * yy(kk) * zz(kk)
     end do

  else if( kfl_paral >= 1 ) then
     !
     ! Parallel: Slaves
     !
     if(nbvar==1) then
        do ii=1,npoi1
           sumxx  = sumxx  + xx(ii) * yy(ii) * zz(ii)
        end do
        do ii=npoi2,npoi3
           sumxx  = sumxx  + xx(ii) * yy(ii) * zz(ii)
        end do
     else
        do ii=1,nbvar*npoi1
           sumxx = sumxx + xx(ii) * yy(ii) * zz(ii)
        end do
        do ii=(npoi2-1)*nbvar+1,npoi3*nbvar
           sumxx  = sumxx  + xx(ii) * yy(ii) * zz(ii)
        end do
     end if
  end if

  if( kfl_paral >= 0 ) then
     !
     ! Parallel: reduce sum
     !
     nparr        =  1
     dummr_par(1) =  sumxx
     parre        => dummr_par
     call par_operat(3_ip)           
     sumxx        =  dummr_par(1)
  end if

end subroutine proxyz

