!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_redieq
  !-----------------------------------------------------------------------
  !****f* Levels/lev_redieq
  ! NAME
  !    lev_redieq
  ! DESCRIPTION
  !    Compute the level set function redistanciation with Sussman equation
  ! USES
  !
  ! USED BY
  !    lev_iniunk
  !    lev_updunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  use def_solver
  use mod_gradie
  use mod_messages, only : livinf
  implicit none
  integer(ip)             :: ipoin

  if( (ittim>=nstre_lev) .and. (nstre_lev/=0) ) return  ! used to stop redistanciation at a certain step to capture possible errors

  call livinf(59_ip,'REDISTANTIATION WITH EQUATION',0_ip)
  if( INOTSLAVE ) print*,' tyred_lev ',tyred_lev

  if ( tyred_lev == 2 ) then
     call lev_redic2()   ! redistance cut nodes only, also returns icupt_lev, searching surfaces only in the element to which the node belongs
  else if ( tyred_lev == 3 .or. tyred_lev == 5 .or. tyred_lev == 6 ) then
!hhelseap     call lev_redisc()   ! redistance cut nodes only, also returns icupt_lev, searching surfaces only in all elements
  else if ( tyred_lev == 4 ) then
     call runend('LEV_REDISC: lev_rediha eliminated in alya 523, it was in the file lev_redisc')
  end if
  
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        flev0_lev(ipoin) = fleve(ipoin,1)
     enddo
  end if

  if(INOTMASTER) then
     call lev_updunk(7_ip) ! UNKNO <= FLEVE
  endif

  if ( tyred_lev == 2 .or. tyred_lev == 3 ) then
     call lev_redsus
  elseif ( tyred_lev == 5 ) then
     call lev_distuc
  elseif ( tyred_lev == 6 ) then
     call lev_distuc
     call lev_redsus
  end if

end subroutine lev_redieq
