!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_roughn()
  !-----------------------------------------------------------------------
  !****f* Domain/roughn
  ! NAME
  !    roughn
  ! DESCRIPTION
  !    Point roughness
  ! OUTPUT
  !    ROUGH
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use mod_memchk
  use def_domain
  use mod_messages, only : messages_live
  use mod_communications_global, only : PAR_SUM
  implicit none
  integer(ip) :: ierro,ipoin,dummi(2)
  real(rp)    :: dummr(2),xdire(3)

  if( kfl_rough /= -1 ) then

     ierro = 0

     if( INOTMASTER ) then

        if( kfl_rough == 0 ) then 
           !
           ! Constant roughness
           !
           do ipoin = 1,npoin
              rough(ipoin) = rough_dom
           end do

        else if( kfl_rough > 0 ) then
           !
           ! Given by a field
           !
           if( nfiel < kfl_rough ) then
              ierro =  1
           else
              rough => xfiel(kfl_rough) % a(1,:,1)
           end if

        end if

     end if
     !
     ! Check errors
     !
     call PAR_SUM(ierro)
     if( ierro /= 0 ) then
        call runend('FILED FOR ROUGHNESS HAS NOT BEEN DEFINED')
     end if
     !
     ! Roughness extension
     !
     if( kfl_rough /= 0 .and. kfl_extro > 0 ) then

        call messages_live('KERMOD: EXTEND ROUGHNESS FROM THE WALL')        
        solve_sol => solve(1:) 
        call inisol()

        if( IMASTER ) then
           call extens(1_ip,dummi,dummr,dummr)
        else  
           xdire        = 0.0_rp
           xdire(ndime) = 1000000.0_rp
           call extens(1_ip,kfl_fixno_rough_ker,xdire,rough)
        end if

        if( ndivi > 0 ) then
           call runend('ROUGHN: CANNOT BE USED BY DIVISOR BECAUSE ROUGHN NOT ZERO IN INTERIOR')
        end if

     end if

  end if

end subroutine ker_roughn

