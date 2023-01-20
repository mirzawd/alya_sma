!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine prodts(nbvar,nbnodes,tt,ss,sumtt,sumts)
  !------------------------------------------------------------------------
  !****f* solite/prodts
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
  use def_kintyp,         only :  ip,rp
  use def_master,         only :  kfl_paral,npoi3
  use mod_communications, only : PAR_SUM
  use def_solver,         only : solve_sol        
  implicit none
  integer(ip), intent(in)  :: nbvar,nbnodes
  real(rp),    intent(in)  :: tt(*),ss(*)
  real(rp),    intent(out) :: sumtt,sumts
  real(rp)                 :: dummr_par(2) 
  integer(ip)              :: kk,ii,jj

  sumtt = 0.0_rp
  sumts = 0.0_rp

  if( solve_sol(1) % kfl_mask == 0 ) then
     
     if(kfl_paral==-1) then
        !
        ! Sequential
        !
        do kk=1,nbvar*nbnodes
           sumtt = sumtt + tt(kk) * tt(kk)
           sumts = sumts + tt(kk) * ss(kk)
        end do
        
     else if(kfl_paral>=1) then
        !
        ! Parallel: Slaves
        !
        sumtt = dot_product(tt(1:nbvar*npoi3),tt(1:nbvar*npoi3))
        sumts = dot_product(tt(1:nbvar*npoi3),ss(1:nbvar*npoi3))
     end if
     
  else

     sumtt = 0.0_rp
     kk = 0
     do ii = 1,npoi3
        do jj = 1,nbvar
           kk = kk + 1
           sumtt = sumtt + tt(kk)*tt(kk)*solve_sol(1) % mask(ii)
           sumts = sumts + tt(kk)*ss(kk)*solve_sol(1) % mask(ii)
        end do
     end do

  end if
  
  if(kfl_paral>=0) then
     !
     ! Parallel: reduce sum
     !
     dummr_par(1) =  sumtt
     dummr_par(2) =  sumts
     call PAR_SUM(2_ip,dummr_par)
     sumtt        =  dummr_par(1)
     sumts        =  dummr_par(2)
  end if

end subroutine prodts

