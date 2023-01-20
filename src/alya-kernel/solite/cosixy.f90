!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine cosixy(nbvar,nbnodes,tt,ss,costs)
  !------------------------------------------------------------------------
  !****f* solite/cosixy
  ! NAME 
  !    prodxy
  ! DESCRIPTION
  !    This routine computes the sclara product of two vectors and the
  !    norm of one of these:
  !    SUMTT = TT.TT
  !    SUMTS = TT.SS
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr
  implicit none
  integer(ip), intent(in)  :: nbvar,nbnodes
  real(rp),    intent(in)  :: tt(*),ss(*)
  real(rp),    intent(out) :: costs
  real(rp)                 :: sumts,sumtt,sumss
  real(rp),    target      :: dummr_par(3) 
  integer(ip)              :: kk,ii,ll

  sumtt = 0.0_rp
  sumss = 0.0_rp
  sumts = 0.0_rp

  if(kfl_paral==-1) then
     !
     ! Sequential
     !
     do kk=1,nbvar*nbnodes
        sumtt = sumtt + tt(kk) * tt(kk)
        sumts = sumts + tt(kk) * ss(kk)
        sumss = sumss + ss(kk) * ss(kk)
     end do

  else if(kfl_paral>=1) then
     !
     ! Parallel: Slaves
     !
     if(nbvar==1) then
        do ii=1,npoi1
           sumtt  = sumtt  + tt(ii) * tt(ii)
           sumts  = sumts  + tt(ii) * ss(ii)
           sumss  = sumss  + ss(ii) * ss(ii)
        end do
        do ii=npoi2,npoi3
           sumtt  = sumtt  + tt(ii) * tt(ii)
           sumts  = sumts  + tt(ii) * ss(ii)
           sumss  = sumss  + ss(ii) * ss(ii)
        end do
     else
        ll=0
        do ii=1,npoi1
           do kk=1,nbvar
              ll=ll+1
              sumtt  = sumtt  + tt(ll) * tt(ll)
              sumts  = sumts  + tt(ll) * ss(ll)
              sumss  = sumss  + ss(ll) * ss(ll)
           end do
        end do
        ll=(npoi2-1)*nbvar
        do ii=npoi2,npoi3
           do kk=1,nbvar
              ll=ll+1
              sumtt  = sumtt  + tt(ll) * tt(ll)
              sumts  = sumts  + tt(ll) * ss(ll)
              sumss  = sumss  + ss(ll) * ss(ll)
           end do
        end do
     end if
  end if

  if(kfl_paral>=0) then
     !
     ! Parallel: reduce sum
     !
     nparr        =  3
     dummr_par(1) =  sumtt
     dummr_par(2) =  sumts
     dummr_par(3) =  sumss
     parre        => dummr_par
     call par_operat(3_ip)           
     sumtt        =  dummr_par(1)
     sumts        =  dummr_par(2)
     sumss        =  dummr_par(3)
  end if

  if( sumss /= 0.0_rp .and. sumtt /= 0.0_rp ) then
     costs = sumts / sqrt(sumtt*sumss)
  else
     costs = 0.0_rp
  end if

end subroutine cosixy

