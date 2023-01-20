!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine resvec(nbvar,nbnodes,tt,ss,difts)
  !------------------------------------------------------------------------
  !****f* solite/resvec
  ! NAME 
  !    resvec
  ! DESCRIPTION
  !    This routine computes the L2 norm:
  !    DIFTS = || TT - SS || / || TT ||
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  kfl_paral,npoi1,npoi2,npoi3,parre,nparr
  implicit none
  integer(ip), intent(in)  :: nbvar,nbnodes
  real(rp),    intent(in)  :: tt(*),ss(*)
  real(rp),    intent(out) :: difts
  real(rp),    target      :: dummr_par(2) 
  real(rp)                 :: xdiff,sumtt
  integer(ip)              :: kk,ii,ll

  difts = 0.0_rp
  sumtt = 0.0_rp

  if(kfl_paral==-1) then
     !
     ! Sequential
     !
     do kk=1,nbvar*nbnodes
        xdiff = tt(kk) - ss(kk)
        difts = difts + xdiff  * xdiff
        sumtt = sumtt + tt(kk) * tt(kk)
     end do

  else if(kfl_paral>=1) then
     !
     ! Parallel: Slaves
     !
     if(nbvar==1) then
        do ii=1,npoi1
           xdiff  = tt(ii) - ss(ii)
           difts  = difts  + xdiff  * xdiff
           sumtt  = sumtt  + tt(ii) * tt(ii)
        end do
        do ii=npoi2,npoi3
           xdiff  = tt(ii) - ss(ii)
           difts  = difts  + xdiff  * xdiff
           sumtt  = sumtt  + tt(ii) * tt(ii)
        end do
     else
        ll=0
        do ii=1,npoi1
           do kk=1,nbvar
              ll=ll+1
              xdiff  = tt(ll) - ss(ll)
              difts  = difts  + xdiff  * xdiff
              sumtt  = sumtt  + tt(ll) * tt(ll)
           end do
        end do
        ll=(npoi2-1)*nbvar
        do ii=npoi2,npoi3
           do kk=1,nbvar
              ll=ll+1
              xdiff  = tt(ll) - ss(ll)
              difts  = difts  + xdiff  * xdiff
              sumtt  = sumtt  + tt(ll) * tt(ll)
           end do
        end do
     end if
  end if

  if(kfl_paral>=0) then
     !
     ! Parallel: reduce sum
     !
     nparr        =  2
     dummr_par(1) =  difts
     dummr_par(2) =  sumtt
     parre        => dummr_par
     call par_operat(3_ip)           
     difts        =  dummr_par(1)
     sumtt        =  dummr_par(2)
  end if

  if(sumtt>0.0_rp) difts=sqrt(difts/sumtt)

end subroutine resvec

