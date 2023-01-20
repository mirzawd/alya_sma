!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_explic
  !-----------------------------------------------------------------------
  !****f* Levels/lev_explic
  ! NAME 
  !    lev_explic
  ! DESCRIPTION
  !    Explicit advance
  ! USES
  ! USED BY
  !    lev_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_levels
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: dtn
  !
  ! Lu^{n-1}
  !
  call bcsrax(1_ip,npoin,1_ip,amatr,c_sol,r_sol,fleve(1,3),unkno)
  !
  ! Service Parall
  !
  call lev_parall(3_ip) ! Lu^{n-1}
  call lev_parall(4_ip) ! RHS

  dtn   =  1.0_rp/dtinv_lev 

  if(kfl_tisch_lev==4) then  


     if(kfl_tiacc_lev==1) then
        !
        ! Runge Kutta 1 
        !
        do ipoin=1,npoin
           unkno(ipoin)=fleve(ipoin,3)+rhsid(ipoin)/vmass(ipoin)
        end do

     end if

     !
     ! BDF2
     !
     !    call lev_elmope(1_ip)


     !    do ipoin=1,npoin
     !       unkno(ipoin)=(-pabdf_lev(2)*fleve(ipoin,3)-pabdf_lev(3)*fleve(ipoin,4)&
     !            +rhsid(ipoin)/vmass(ipoin))/pabdf_lev(1)
     !    end do

     !  end if

  end if
  !
  ! Impose boundary conditions
  !
  if(kfl_advec_lev==4) then
     if(kfl_timet_lev==1) then
        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)>0) unkno(ipoin)=fleve(ipoin,3)-dtn
           !   if(kfl_fixno_lev(ipoin)>0) unkno(ipoin)=bvess_lev(ipoin)
        end do
     end if
  else if(kfl_advec_lev==6) then
     if(kfl_timet_lev==1) then
        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)>0) unkno(ipoin)=fleve(ipoin,3)
        end do
     end if
  end if

end subroutine lev_explic
