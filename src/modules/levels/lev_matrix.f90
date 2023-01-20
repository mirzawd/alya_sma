!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_matrix
  !-----------------------------------------------------------------------
  !****f* Levels/lev_matrix
  ! NAME 
  !    lev_matrix
  ! DESCRIPTION
  !    Compute matrix and RHS
  ! USES
  !    lev_elmope
  ! USED BY
  !    lev_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  use def_solver
  implicit none
  integer(ip)       :: izmat,ipoin

  if(kfl_timet_lev==1) then ! Explicit treatment case 

     do ipoin=1,npoin
        rhsid(ipoin)=0.0_rp
     end do
     do izmat=1,solve(1)%nzmat
        amatr(izmat)=0.0_rp
     end do
     !
     ! LHS assembly
     !     
     call lev_elmope(1_ip)

  else if(kfl_timet_lev==2) then ! Implicit treatment case 

     do izmat=1,solve(1)%nzmat
        amatr(izmat)=0.0_rp
     end do
     do ipoin=1,solve(1)%nzrhs
        rhsid(ipoin)=0.0_rp
     end do
     call lev_elmope(2_ip)
  endif
  !
  ! Penalization
  !
  !thick_lev = 0.06_rp
  !do ipoin = 1,npoin
  !   izdod = r_dom(ipoin) - 1
  !   jpoin = 0
  !   do while( jpoin /= ipoin )
  !      izdod = izdod + 1
  !      jpoin = c_dom(izdod)
  !   end do
  !   lambda = 0.9_rp/thick_lev
  !   eps    = 0.5_rp-0.5_rp*tanh(lambda*dista_lev(ipoin))**2
  !   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
  !      amatr(izdom) = amatr(izdom) * (1.0_rp-eps)
  !   end do
  !   rhsid(ipoin) = rhsid(ipoin) * (1.0_rp-eps) 
  !   if( kfl_fixno_lev(1,ipoin) == 0 ) then
  !      amatr(izdod) = amatr(izdod) + eps * vmass(ipoin)
  !      rhsid(ipoin) = rhsid(ipoin) + eps * vmass(ipoin) * dista_lev(ipoin)
  !   end if
  !end do

end subroutine lev_matrix
