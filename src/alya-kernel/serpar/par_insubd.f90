!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_insubd(ii,jpoin,ifoun)
  !------------------------------------------------------------------------
  !****f* parall/par_insubd
  ! NAME 
  !    par_insubd
  ! DESCRIPTION
  !    Check if point JPOIN is in boundary with my neighbor II
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip
  use mod_parall, only     :  commd
  implicit none
  integer(ip), intent(in)  :: ii,jpoin
  integer(ip), intent(out) :: ifoun
  integer(ip)              :: ll

  ifoun = 0
  ll    = commd % bound_size(ii) 
  do while( ll <= commd % bound_size(ii+1)-1 )
     if( commd % bound_perm(ll) == jpoin ) then
        ifoun = 1
        return
     end if
     ll = ll + 1
  end do
  
end subroutine par_insubd
